#include <pybindings.h>
#include <container_pybindings.h>

#include <frbutils/frbutils.h>
#include <G3Logging.h>
#include <G3Vector.h>
#include <G3Data.h>
#include <G3Map.h>
#include <G3Vector.h>

#include <G3Logging.h>

#include <iostream>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <todfilter/numericutils.h>
#include <map>
#include <set>
#include <stack>
#include <unordered_set>
#include <algorithm>

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif


G3_SERIALIZABLE_CODE(FrbDetInfo);
G3_SERIALIZABLE_CODE(G3VectorFrbDetInfo);

G3_SERIALIZABLE_CODE(FrbEventInfo);
G3_SERIALIZABLE_CODE(G3VectorFrbEventInfo);

double get_poisson_threshold( double total_num_counts,
			      double num_types,
			      double probability_cutoff){
	if (num_types == 0 ) return 0;
	double mean = total_num_counts/num_types;
	boost::math::poisson_distribution<double> pd(mean);
	return boost::math::quantile(pd, probability_cutoff);
}



double get_log_likelihood_from_model_svd( const GslVectorPtr y_in, double variance, 
					  const GslMatrixPtr mod,
					  const GslMatrixPtr U, const GslMatrixPtr V, 
					  const GslVectorPtr S, 
					  std::vector<double> & out_fit_params){
	g3_assert(variance > 0);
	double output = 0.0;  
	size_t M = U->size1;
	size_t N = U->size2;
	gsl_vector * y = y_in.get();
	gsl_vector * x = gsl_vector_alloc(N);
	gsl_vector * y_fit = gsl_vector_alloc(M);
	int status = gsl_linalg_SV_solve (U.get(), V.get(), S.get(), y, x);
	gsl_blas_dgemv (CblasNoTrans, 1.0, mod.get(), x, 0, y_fit);
	for (size_t i=0; i < M; i++){
		double err = y->data[i] - y_fit->data[i];
		output += -0.5 * err * err / variance;
	}
  
	g3_assert(out_fit_params.size() == N);
	memcpy(&(out_fit_params[0]), x->data, N * sizeof(double));
	
	gsl_vector_free(x);
	gsl_vector_free(y_fit);
	return output;
}



GslMatrixPtr construct_poly_delta_heaviside_model(size_t m_len, size_t poly_order, 
						  bool include_heavi, bool include_delta){
	log_debug("Enter");
	size_t n_params = poly_order + include_heavi + include_delta + 1;
	GslMatrixPtr m = get_shared_gsl_matrix( m_len, n_params);
	double * out_model_matrix = m->data;
	//poly or hpf always have a mean filter
	for (size_t i=0; i < m_len; i++){
		out_model_matrix[n_params * i] = 1.0f;
	} 
	//cover the first case
	if (poly_order >0){
		int j = 1;
		for (size_t i=0; i < m_len; i++){
			double x = -1.0f + 2 * ((double)i)/((double) m_len - 1.0f);
			out_model_matrix[n_params * i + j] = x;
		}
	}
	for (size_t i=0; i < m_len; i++){
		for (size_t j = 2; j <  1 + poly_order; j++){
			double x = -1.0f + 2 * ((double)i)/((double) m_len - 1.0f);
			double j_leg = j;
			out_model_matrix[n_params * i + j] = ((2*j_leg-1) * x * out_model_matrix[(j-1)+i*n_params] -
							      (j_leg-1) * out_model_matrix[(j-2)+i*n_params])/j_leg;
		}
	}
	
	size_t offset = poly_order+1;
	
	if (include_heavi){
		size_t heav_loc = m_len/2;
		for (size_t i=0; i < heav_loc; i++) out_model_matrix[n_params * i + offset] = -0.5;
		for (size_t i=heav_loc+1; i < m_len; i++) out_model_matrix[n_params * i + offset] = 0.5;
		offset++;
	}
	
	if (include_delta){
		size_t delt_loc = m_len/2;
		out_model_matrix[n_params * delt_loc + offset] = 1.0;
		offset++;
	}
	log_debug("Exit");
	return m;
}


void get_ts_poly_delta_heavi_ll(G3TimestreamMapConstPtr ts_map_const, 
				G3MapDoubleConstPtr variance,
				size_t fit_model_len, size_t poly_order, 
				bool include_heavi, bool include_delta,
				G3TimestreamMapPtr ll_map, 
				G3TimestreamMapPtr amp_map,
				G3TimestreamMapPtr heavi_amp_map,
				int special_index = -1
	){
	log_debug("Enter");
	//This is actually const, just gsl, mreh
	G3TimestreamMap * ts_map = const_cast<G3TimestreamMap *>(ts_map_const.get());
	//makes certain model is not too large
	size_t ts_len = ts_map->begin()->second->size();
	g3_assert(ts_len >= fit_model_len);
	//constructs model svd  
	GslMatrixPtr mod = construct_poly_delta_heaviside_model(fit_model_len, poly_order, 
								include_heavi, include_delta);
	
	log_debug("Pdelta\n");
	const size_t start_offset = fit_model_len/2;
	const size_t stop_offset = fit_model_len/2 + fit_model_len%2;
	
	size_t ts_iteration_start = start_offset;
	size_t ts_iteration_stop = ts_len - stop_offset;
	
	if (special_index >= 0 ) {
		g3_assert(special_index >= ts_iteration_start &&
			  special_index < ts_iteration_stop);
		ts_len = 1;
		ts_iteration_start = special_index;
		ts_iteration_stop = special_index + 1;
	}

	log_debug("Thin svd\n");
	GslMatrixPtr U,V;
	GslVectorPtr S;
	get_thin_svd_of_matrix(mod, U, V, S);
	
	std::vector<std::string> id_lst(ts_map->size());
	size_t j = 0;
	for (auto it = ts_map->begin(); it != ts_map->end(); it++){
		(*ll_map)[it->first] = G3TimestreamPtr(new G3Timestream(ts_len, 0));
		if (include_delta) 
			(*amp_map)[it->first] = G3TimestreamPtr(new G3Timestream(ts_len, 0));
		if (include_delta && include_heavi) 
			(*heavi_amp_map)[it->first] =G3TimestreamPtr(new G3Timestream(ts_len, 0));
		id_lst[j] = it->first;
		j++;
	}
	
#ifdef OPENMP_FOUND
#pragma omp parallel for
#endif
	for (j=0; j < id_lst.size(); j++){
		std::vector<double> fit_params(S->size);
		GslVectorPtr ts_data = GslVectorPtr(new gsl_vector);
		ts_data->size = fit_model_len;
		ts_data->stride = 1;
		ts_data->block = NULL;
		ts_data->owner = 0;
		G3TimestreamPtr ll = ll_map->at(id_lst[j]);
		G3TimestreamPtr amp;
		G3TimestreamPtr hamp;
		if (include_delta) amp = G3TimestreamPtr(new G3Timestream(ts_len, 0));
		if (include_delta && include_heavi) hamp = G3TimestreamPtr(new G3Timestream(ts_len, 0));
		for (size_t i = ts_iteration_start; i < ts_iteration_stop; i++){
			size_t store_ind = i;
			if (special_index >= 0) {
				store_ind = 0;
			} 

			ts_data->data = &( (*ts_map->at(id_lst[j]))[i-start_offset]);
			(*ll)[store_ind] = get_log_likelihood_from_model_svd( 
				ts_data, variance->at(id_lst[j]), mod, U, V, S, fit_params);
			if (include_delta) 
				(*amp)[store_ind] = fit_params[mod->size2-1];
			if (include_delta && include_heavi) 
				(*hamp)[store_ind] = fit_params[mod->size2-2];
		}

		if (include_delta)   
			(*amp_map)[id_lst[j]] = amp;
		if (include_delta && include_heavi)  
			(*heavi_amp_map)[id_lst[j]] = hamp;
	}
	
	log_debug("Exit");
}


void search_for_trigger( G3TimestreamMapConstPtr ll_model,  G3TimestreamMapConstPtr ll_baseline,
			 double trigger_threshold, size_t min_distance,
			 std::vector< size_t > & triggered_events){
	log_debug("search_for_trigger");
	triggered_events.clear();
	std::unordered_set<size_t> full_hits;
	for (auto it = ll_model->begin(); it != ll_model->end(); it++){
		const G3Timestream & ll_m = *ll_model->at(it->first);
		const G3Timestream & ll_b = *ll_baseline->at(it->first);
		G3Timestream ll_delta = ll_m-ll_b;
		for (size_t i=0; i < ll_m.size(); i++){
			if (ll_delta[i] > trigger_threshold){
				int local_max;
				double max_val;
				find_local_maxima<double>(ll_delta, 
							  i - min_distance, 
							  i + min_distance,
							  local_max, max_val);
				full_hits.insert(local_max);
			}
		}
	}
	triggered_events = std::vector<size_t> (full_hits.begin(), full_hits.end());
	log_debug("Exit");
}


void look_for_other_bolos_seeing_it(const std::vector< size_t > & triggered_events_start, 
				    G3TimestreamMapConstPtr ll_model,  
				    G3TimestreamMapConstPtr ll_baseline,
				    double alt_trigger_threshold,
				    G3TimestreamMapConstPtr fit_amp,
				    G3TimestreamMapConstPtr fit_hamp,
				    G3VectorFrbEventInfoPtr events_ptr, size_t search_width){
	log_debug("Enter");
	G3VectorFrbEventInfo & events = *events_ptr;
	
	std::vector< size_t > triggered_events(triggered_events_start.begin(), 
					       triggered_events_start.end());
	std::sort(  triggered_events.begin(), triggered_events.end());
	
	
	std::vector< size_t > search_starts(triggered_events.size(), 0);
	std::vector< size_t > search_stops(triggered_events.size(), 0);
	std::vector< bool > ignore_ev(triggered_events.size(), false);
	
	for (size_t j = 0; j < triggered_events.size(); ){
		search_starts[j] = triggered_events[j] < search_width ? 
							 0 : triggered_events[j] - search_width;
		search_stops[j]  = triggered_events[j] + search_width;
		search_stops[j] = search_stops[j] > ll_model->begin()->second->size() - 1 ?
			ll_model->begin()->second->size() - 1 : search_stops[j];
		ignore_ev[j] = false;
		
		// look to see if our search width in encompassing another event.  
		//If so we skip that event and increase this ones search size.
		size_t tmp_i = j+1;
		while ( tmp_i < triggered_events.size() ){
			ignore_ev[tmp_i] = true;
			if (triggered_events[tmp_i] <= search_stops[j]){
				search_stops[j] = triggered_events[tmp_i] + search_width;
				search_stops[j] = search_stops[j] > ll_model->begin()->second->size() - 1 ?
					ll_model->begin()->second->size() - 1 : search_stops[j];
			} else {
				break;
			}
			tmp_i++;
		}
		j = tmp_i;
	}      
	
	for (size_t j=0; j < triggered_events.size(); j++) 
		g3_assert(search_stops[j] < ll_model->begin()->second->size() );
	
	for (size_t j = 0; j < triggered_events.size(); j++){
		if ( ignore_ev[j] ) continue;
		
		size_t * it = &(triggered_events[j]);
		FrbEventInfo ev_info;
		
		ev_info.scan_index = *it;

		size_t min_ind = search_starts[j];
		size_t max_ind = search_stops[j];
		
		for (auto det_it = ll_model->begin(); det_it != ll_model->end(); det_it++){
			int largest_ind = *it;
			double largest_val = (*ll_model->at(det_it->first))[*it] - 
				(*ll_baseline->at(det_it->first))[*it];
			
			//g3_assert(largest_val >= 0);

			//search through all samples in a range
			for (size_t i = min_ind; i <= max_ind; i++){
				double llm = (*ll_model->at(det_it->first))[i];
				double llb = (*ll_baseline->at(det_it->first))[i];
				if ((llm - llb) > largest_val) {
					largest_val = llm - llb;
					largest_ind = i;
				}
			}
			//if one of the samples is large enough, use that
			if (largest_val > alt_trigger_threshold){
				FrbDetInfo det_info;
				det_info.bid = det_it->first;

				det_info.ll_model = (*ll_model->at(det_it->first))[largest_ind];
				det_info.ll_baseline = (*ll_baseline->at(det_it->first))[largest_ind];
				det_info.significance = det_info.ll_model - det_info.ll_baseline;
				det_info.amplitude = (*fit_amp->at(det_it->first))[largest_ind];
				det_info.heavi_amp = (*fit_hamp->at(det_it->first))[largest_ind];


				det_info.scan_index = largest_ind;
				ev_info.det_info.push_back(det_info);
			}
		}
		events.push_back(ev_info);
	}
	log_debug("Exit");
}




PolyLikelihoodFiller::PolyLikelihoodFiller(size_t model_len, size_t poly_order, 
					   bool include_heaviside, 
					   bool include_delta, 
					   std::string ts_key, 
					   std::string variance_key,
					   std::string loglike_output_key,
					   std::string amp_map_output_key,
					   std::string hamp_map_output_key
	) :
	model_len_(model_len), poly_order_(poly_order), include_heaviside_(include_heaviside),
	include_delta_(include_delta), ts_key_(ts_key), variance_key_(variance_key),
	loglike_output_key_(loglike_output_key), amp_map_output_key_(amp_map_output_key),
	hamp_map_output_key_(hamp_map_output_key)
{
	g3_assert(model_len % 2 == 1);
}

void PolyLikelihoodFiller::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	log_debug("Enter");
	if (frame->type == G3Frame::Scan){
		G3TimestreamMapConstPtr ts_map = frame->Get<G3TimestreamMap>(ts_key_);
		G3MapDoubleConstPtr variance = frame->Get<G3MapDouble>(variance_key_);
    
		G3TimestreamMapPtr ll_map(new G3TimestreamMap);
		G3TimestreamMapPtr amp_map(new G3TimestreamMap);
		G3TimestreamMapPtr hamp_map(new G3TimestreamMap);
		
		get_ts_poly_delta_heavi_ll(ts_map, variance,
					   model_len_, poly_order_, 
					   include_heaviside_, include_delta_,
					   ll_map, 
					   amp_map, hamp_map);
		frame->Put(loglike_output_key_, ll_map);
		if (include_delta_) frame->Put(amp_map_output_key_, amp_map);
		if (include_delta_ && include_heaviside_) frame->Put(hamp_map_output_key_, hamp_map);
	} 
	out.push_back(frame);
	log_debug("Exit");
}

DeltaEventHunter::DeltaEventHunter(std::string ll_model_key,  
				   std::string ll_base_key,  
				   double trigger_thresh, 
				   double other_det_thresh,
				   size_t min_distance,
				   std::string output_event_key,
				   std::string fit_amp_key,
				   std::string fit_hamp_key,
				   int search_width
	):
	ll_model_key_(ll_model_key), ll_base_key_(ll_base_key),
	trigger_thresh_(trigger_thresh),
	other_det_thresh_(other_det_thresh), min_distance_(min_distance),
	output_event_key_(output_event_key), fit_amp_key_(fit_amp_key), fit_hamp_key_(fit_hamp_key),
	search_width_(search_width)
{}

void DeltaEventHunter::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	log_debug("Enter");
	if (frame->type == G3Frame::Scan){
		G3TimestreamMapConstPtr ll_mod = frame->Get<G3TimestreamMap>(ll_model_key_);
		G3TimestreamMapConstPtr ll_bas = frame->Get<G3TimestreamMap>(ll_base_key_);
		G3TimestreamMapConstPtr fit_amp = frame->Get<G3TimestreamMap>(fit_amp_key_);
		G3TimestreamMapConstPtr fit_hamp = frame->Get<G3TimestreamMap>(fit_hamp_key_);
		
		std::vector<size_t> triggered_events;
		G3VectorFrbEventInfoPtr events(new G3VectorFrbEventInfo);
		search_for_trigger( ll_mod, ll_bas,  trigger_thresh_, 
				    min_distance_, triggered_events);
		look_for_other_bolos_seeing_it(triggered_events, ll_mod, ll_bas, 
					       other_det_thresh_, fit_amp, fit_hamp,
					       events, search_width_);
		frame->Put(output_event_key_, events);
	} 
	out.push_back(frame);
	log_debug("Exit");
}




int get_hist_bin(const std::string & main_det, 
		 const std::string & alt_det,
		 const G3MapDouble & tes_x_pos,
		 const G3MapDouble & tes_y_pos,
		 const G3VectorDouble & histogram_bins) {

	int hist_bin;
	double main_x = tes_x_pos.at(main_det);
	double main_y = tes_y_pos.at(main_det);
	double alt_x = tes_x_pos.at(alt_det);
	double alt_y = tes_y_pos.at(alt_det);

	double dist = ( (main_x-alt_x)*(main_x-alt_x) + 
			(main_y-alt_y)*(main_y-alt_y)   );
	hist_bin = histogram_bins.size() -1;
	for (size_t i = 0; i < histogram_bins.size() - 1; i++){
	  if (dist >= ( histogram_bins.at(i) * histogram_bins.at(i) )  &&
	      dist < ( histogram_bins.at(i+1) * histogram_bins.at(i+1) ) ) {
			hist_bin = i;
			break;
		}
	}

	return hist_bin;
}

void get_exposure_hist(const std::string & main_det, 
		       const std::vector<std::string> & det_ids,
		       const G3MapDouble & tes_x_pos,
		       const G3MapDouble & tes_y_pos,
		       const G3VectorDouble & histogram_bins,
		       G3VectorDouble & exposure_hist
	) {
	G3VectorDouble bins(histogram_bins);
	exposure_hist = G3VectorDouble(histogram_bins.size(), 0);

	//printf("Bins size %zu\nExpose %zu\n", bins.size(), exposure_hist.size());

	double main_x = tes_x_pos.at(main_det);
	double main_y = tes_y_pos.at(main_det);

	for (size_t i=0; i < bins.size(); i++) bins[i] = bins[i] * bins[i];

	for (auto it = det_ids.begin(); it != det_ids.end(); it++) {
		double delta_x = main_x - tes_x_pos.at(*it);
		double delta_y = main_y - tes_y_pos.at(*it);
		double dist2 = delta_x * delta_x + delta_y * delta_y;
		size_t bin = bins.size() - 1;
		for (size_t i=0; i < bins.size() - 1; i++) {
			if ( dist2 >= bins[i] && dist2 < bins[i+1]) {
				bin = i;
				break;
			}
		}
		exposure_hist[bin] += 1;
	}
}


void get_exposure_hist_even_bin(const std::string & main_det, 
				const std::vector<std::string> & det_ids,
				const G3MapDouble & tes_x_pos,
				const G3MapDouble & tes_y_pos,
				const G3VectorDouble & histogram_bins,
				G3VectorDouble & exposure_hist
	) {

  //maybe works?
	G3VectorDouble bins(histogram_bins);
	exposure_hist = G3VectorDouble(histogram_bins.size(), 0);


	double main_x = tes_x_pos.at(main_det);
	double main_y = tes_y_pos.at(main_det);

	double min = histogram_bins[0];
	double max = histogram_bins[bins.size()-1];
	double delt = max - min;
	// (val-min) /( max - min)
	for (auto it = det_ids.begin(); it != det_ids.end(); it++) {
		double delta_x = main_x - tes_x_pos.at(*it);
		double delta_y = main_y - tes_y_pos.at(*it);
		double dist2 = sqrt(delta_x * delta_x + delta_y * delta_y);
		size_t bin;
		if ( dist2 < min || dist2 > max) {
			bin = histogram_bins.size()-1;
		} else {
		  bin = (dist2 - min) / delt * (histogram_bins.size() -1);
		}
		exposure_hist[bin] += 1;
	}
}


G3VectorFrbEventInfo filter_frb_events(
	const G3VectorFrbEventInfo in_ev,
	const std::vector<int> is_bad) {
	G3VectorFrbEventInfo out_ev(in_ev);

	if (in_ev.size() == 0) return out_ev;
	g3_assert(is_bad.size() == out_ev.size());
	size_t n_bad = 0;
	size_t tsize = out_ev.size();
	for (int i = is_bad.size() - 1; i >= 0; i--) {
		if (is_bad[i]) {
			n_bad++;
			/**
			printf("%zu %zu %zu\n", i, out_ev.size(), is_bad.size());
			g3_assert(tsize - n_bad >= 0);
			g3_assert(i >= 0);
			g3_assert(tsize - n_bad < out_ev.size());
			g3_assert(i <out_ev.size() );
			**/
			out_ev[i] = out_ev[tsize - n_bad];
		}
	}
	out_ev.resize(out_ev.size() - n_bad);
	return out_ev;
}





int get_num_dets_dist_away(const std::string & main_det, 
			   const std::vector<std::string> & det_ids,
			   const G3MapDouble & tes_x_pos,
			   const G3MapDouble & tes_y_pos,
			   const G3VectorDouble & histogram_bins,
			   int hist_bin
	){
	int n_dets;
	double main_x = tes_x_pos.at(main_det);
	double main_y = tes_y_pos.at(main_det);	
	double low_val = histogram_bins.at(hist_bin) * histogram_bins.at(hist_bin);
	double high_val = histogram_bins.at(hist_bin + 1) * histogram_bins.at(hist_bin + 1);
	n_dets = 0;
	for (auto it = det_ids.begin(); it != det_ids.end(); it++){
		double x_delt = main_x - tes_x_pos.at(*it);
		double y_delt = main_y - tes_y_pos.at(*it);
		double dist = x_delt*x_delt + y_delt*y_delt;
		n_dets += ( (dist >= low_val) && 
			    (dist < high_val) );
	}
	return n_dets;
}

int get_num_det_pairs(const G3MapVectorString & det_group_ids_base,
		      const G3MapDouble & tes_x_pos,
		      const G3MapDouble & tes_y_pos,
		      const G3VectorDouble & histogram_bins,
		      const G3VectorString & good_detectors, //consider swapping with unordered_map for o(1) lookup
		      int hist_bin
	){
	double low_val = histogram_bins.at(hist_bin) * histogram_bins.at(hist_bin);
	double high_val = histogram_bins.at(hist_bin + 1) * histogram_bins.at(hist_bin + 1);
	int n_pairs = 0;

	G3MapVectorString det_group_ids(det_group_ids_base);
	//std::vector<std::string> keys;
	//for (auto it = det_group_ids.begin(); it != det_group_ids.end(); it++) keys.push_back(it->first);

	std::set<std::string> good_ds(good_detectors.begin(), good_detectors.end());
	//Removes the bad detectors
	for (auto git = det_group_ids.begin(); git != det_group_ids.end(); git++){
		std::vector<std::string> &  gitr = git->second;
		//for (auto it = keys.begin(); it != keys.end(); it++ ){
		//auto gitr = det_group_ids[*it];
		int n_bad = 0;
		int end_condition = gitr.size();
		for (size_t i = 0; i < end_condition; ){
			//if(std::find(good_detectors.begin(), good_detectors.end(), gitr[i]) != good_detectors.end()) {
			if (good_ds.count(gitr[i])) {
				i++;
			} else {
				gitr[i] = gitr[end_condition-1];
				end_condition--;
				n_bad++;
			}
		}
		for (size_t i=0; i < n_bad; i++){ 
			gitr.pop_back();
		}
	}

	//calculates all the pairs of detectors (including the double counting needed)
	for (auto git = det_group_ids.begin(); git != det_group_ids.end(); git++){
		std::vector<std::string> &  gitr = git->second;
		for (size_t i = 0; i < gitr.size(); i++){
			//log_error("before the strom %s",gitr[i].c_str());
			if (tes_x_pos.count(gitr[i] ) == 0) continue;
			float x = tes_x_pos.at(gitr[i] );
			float y = tes_y_pos.at(gitr[i] );
			for (size_t j=0; j < gitr.size(); j++){
				if (i==j) continue;
				if (tes_x_pos.count(gitr[j] ) == 0) continue;
				double x_delt = x - tes_x_pos.at( gitr[j] );
				double y_delt = y - tes_y_pos.at( gitr[j] );
				double dist = x_delt*x_delt + y_delt*y_delt;
				n_pairs += ( (dist >= low_val) && 
					     (dist < high_val) );
			}
		}
	}
	return n_pairs;
}

namespace bp = boost::python;

PYBINDINGS("frbutils") {
  bp::class_<PolyLikelihoodFiller, bp::bases<G3Module>,
             boost::shared_ptr<PolyLikelihoodFiller>,
             boost::noncopyable >("PolyLikelihoodFiller", bp::init<size_t, size_t,bool, 
				  bool, std::string, std::string, 
				  std::string, std::string, std::string>
                                  ( ( bp::arg("model_len"),
				      bp::arg("poly_order"),
				      bp::arg("include_heaviside"),
				      bp::arg("include_delta"),
				      bp::arg("ts_key"),
				      bp::arg("variance_key"),
				      bp::arg("loglike_output_key"),
				      bp::arg("amp_map_output_key"),
				      bp::arg("hamp_map_output_key")
				      ))
                                  )
    .def_readonly("__g3module__", true)
    .def("Process", &PolyLikelihoodFiller::Process)
    ;
  bp::implicitly_convertible< boost::shared_ptr<PolyLikelihoodFiller>, G3ModulePtr>();



  EXPORT_FRAMEOBJECT(FrbDetInfo, init<>(), "FRB Detector Info")
	  .def_readwrite("bid", &FrbDetInfo::bid)
	  .def_readwrite("significance", &FrbDetInfo::significance)
	  .def_readwrite("ll_model", &FrbDetInfo::ll_model)
	  .def_readwrite("ll_baseline", &FrbDetInfo::ll_baseline)
	  .def_readwrite("amplitude", &FrbDetInfo::amplitude)
	  .def_readwrite("heavi_amp", &FrbDetInfo::heavi_amp)

	  .def_readwrite("board_id", &FrbDetInfo::board_id)
	  .def_readwrite("module_id", &FrbDetInfo::module_id)

	  .def_readwrite("squid_sig", &FrbDetInfo::squid_sig)
	  .def_readwrite("squid_sig_max", &FrbDetInfo::squid_sig_max)
	  .def_readwrite("wafer_sig", &FrbDetInfo::wafer_sig)

	  .def_readwrite("ra", &FrbDetInfo::ra)
	  .def_readwrite("dec", &FrbDetInfo::dec)

	  .def_readwrite("was_injected", &FrbDetInfo::was_injected)
	  .def_readwrite("scan_index", &FrbDetInfo::scan_index)

	  .def_readwrite("q_significance", &FrbDetInfo::q_significance)
	  .def_readwrite("q_ll_model", &FrbDetInfo::q_ll_model)
	  .def_readwrite("q_ll_baseline", &FrbDetInfo::q_ll_baseline)
	  .def_readwrite("q_amplitude", &FrbDetInfo::q_amplitude)
	  .def_readwrite("q_heavi_amp", &FrbDetInfo::q_heavi_amp)
	  ;

  bp::implicitly_convertible< boost::shared_ptr<FrbDetInfo>, G3FrameObjectPtr>();

  register_pointer_conversions<FrbDetInfo>();
  register_vector_of<FrbDetInfo>("FrbDetInfo");
  register_g3vector<FrbDetInfo>("G3VectorFrbDetInfo", "List of FrbDetInfos");

  EXPORT_FRAMEOBJECT(FrbEventInfo, init<>(), "FRB Event Info")
	  .def_readwrite("det_info", &FrbEventInfo::det_info)
	  .def_readwrite("observation_name", &FrbEventInfo::observation_name)
	  .def_readwrite("observation_number", &FrbEventInfo::observation_number)
	  
	  .def_readwrite("scan_number", &FrbEventInfo::scan_number)
	  .def_readwrite("scan_index", &FrbEventInfo::scan_index)

	  .def_readwrite("event_time", &FrbEventInfo::event_time)

	  .def("filter_det_evs", &FrbEventInfo::filter_det_evs)
	  ;
  bp::implicitly_convertible< boost::shared_ptr<FrbEventInfo>, G3FrameObjectPtr>();



  bp::def("get_num_dets_dist_away", get_num_dets_dist_away);
  bp::def("get_hist_bin", get_hist_bin);
  
  bp::class_<DeltaEventHunter, bp::bases<G3Module>,
             boost::shared_ptr<DeltaEventHunter>,
             boost::noncopyable >("DeltaEventHunter", bp::init<std::string, std::string, 
				  double, double, size_t, 
				  std::string, std::string, 
				  std::string, int>
                                  ( ( bp::arg("ll_model_key"),
				      bp::arg("ll_base_key"),
				      bp::arg("trigger_thresh"),
				      bp::arg("other_det_thresh"),
				      bp::arg("min_distance"),
				      bp::arg("output_event_key"),
				      bp::arg("fit_amp_key"),
				      bp::arg("fit_hamp_key"),
				      bp::arg("search_width")				      
				      ))
                                  )
    .def_readonly("__g3module__", true)
    .def("Process", &DeltaEventHunter::Process)
    ;
  bp::implicitly_convertible< boost::shared_ptr<DeltaEventHunter>, G3ModulePtr>();

  register_pointer_conversions<FrbEventInfo>();
  register_vector_of<FrbEventInfo>("FrbEventInfo");
  register_g3vector<FrbEventInfo>("G3VectorFrbEventInfo", "List of candidate FRB events");

  bp::def( "get_num_det_pairs", get_num_det_pairs );
  bp::def( "get_exposure_hist", get_exposure_hist );
  bp::def( "get_exposure_hist_even_bin", get_exposure_hist_even_bin );

  bp::def( "filter_frb_events", filter_frb_events);

  bp::def( "get_ts_poly_delta_heavi_ll", get_ts_poly_delta_heavi_ll);

  bp::def( "get_poisson_threshold", get_poisson_threshold);

} 
