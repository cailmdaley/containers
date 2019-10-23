#include <pyconfig.h>
#include <frbutils/frbfiltering.h>
#include <todfilter/numericutils.h>

#include <G3Data.h>
#include <G3Logging.h>

#include <iostream>

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif



double get_num_threshold( double total_num_counts, 
			  double num_types,
			  double probability_cutoff){
  if (num_types == 0 ) return 0;
  double mean = total_num_counts/num_types;
  boost::math::poisson_distribution<double> pd(mean);
  return boost::math::quantile(pd, probability_cutoff);
}

void get_too_many_frb_events(const G3VectorFrbEventInfoPtr in_vec, 
			     double probability_cutoff, int minimum_for_cut,
			     std::unordered_set<std::string> & bad_det
			     ){
  size_t total_det_hits = 0;
  G3MapInt det_hits;
  
  //sum all of the obs hits and det hits
  const G3VectorFrbEventInfo & iv = *in_vec;
  for (auto it=iv.begin(); it != iv.end(); it++){
    const FrbEventInfo & info = *it;
    for (auto dit=info.triggers.begin(); dit!=info.triggers.end(); dit++){
      if (det_hits.find(*dit) == det_hits.end()) det_hits[*dit] = 1;
      else det_hits[*dit] += 1;      
      total_det_hits++;
    }
  }
  //find the cutoff
  double det_thresh = get_num_threshold( total_det_hits,  det_hits.size(), probability_cutoff);
  det_thresh = MAX(minimum_for_cut, det_thresh);
  printf("Cut is %lf\n", det_thresh);


  //log_warn("Using detector number threshold of %lf with total hits %zu\n", det_thresh, total_det_hits);

  int total_cut = 0;
  for (auto it = det_hits.begin(); it != det_hits.end(); it++){
    //log_warn("%s %d", it->first.c_str(), it->second);
    if (it->second > det_thresh) {
      bad_det.insert(it->first);
      total_cut++;
    }
  }
}




double get_individual_delta_heavi_poly_fit_ll( const G3Timestream & ts,
					       double variance,
					       size_t index_we_want_delta_at,
					       int poly_order, size_t model_len,
					       bool include_delta, bool include_heavi){
  g3_assert(index_we_want_delta_at  >= model_len/2);

  G3Timestream * ts_deconst = const_cast<G3Timestream * > (&ts);
  GslVectorPtr ts_data = GslVectorPtr(new gsl_vector);
  ts_data->size = model_len;
  ts_data->stride = 1;
  ts_data->block = NULL;
  ts_data->owner = 0;
  ts_data->data = &(  ts_deconst->at( index_we_want_delta_at - model_len/2 ) );

  std::vector<double> fit_params(poly_order + include_heavi + include_delta + 1);
  GslMatrixPtr mod;
  GslMatrixPtr U,V;
  GslVectorPtr S;
  mod = construct_poly_delta_heaviside_model(model_len, poly_order, include_heavi, include_delta);
  get_thin_svd_of_matrix(mod, U, V, S);
  return get_log_likelihood_from_model_svd( ts_data, variance, mod, U, V, S, fit_params);
}





double get_group_sum_loglikelihood(const FrbEventInfo & frb_event,
				   const std::vector<std::string> & det_lst, 
				   const G3MapDouble & variance_map,
				   bool exclude_dets_in_event,
				   size_t fit_model_len, size_t poly_order,
				   bool include_heavi, 
				   double var_low_cutoff, double var_high_cutoff
				   ){
  log_debug("entering");
  log_debug("frb event stored ts size %zu", frb_event.event_ts.size());
  
  std::vector<std::string> dlst;
  if (exclude_dets_in_event){
    for (size_t i=0; i < det_lst.size(); i++){
      bool found_det = false;
      for (size_t j=0; j < frb_event.triggers.size(); j++){
	if (det_lst[i] == frb_event.triggers[j]) {
	  found_det = true;
	  break;
	}
      }
      if (!found_det){
	dlst.push_back(det_lst[i]);
      }
    }
  } else{
    dlst = det_lst;
  }
  if (dlst.size() == 0){
    log_error("no detectors surivive unique cut wtf");
    return -1;
  }
  size_t ts_len = frb_event.event_ts.begin()->second->size();
  size_t ev_index = ts_len/2;
  G3Timestream sum_ts(ts_len, 0.0);
  double sum_variance = 0.0;
  for (auto ts = dlst.begin(); ts != dlst.end(); ts++){
    double var =  variance_map.at(*ts);
    if ( variance_map.find(*ts) !=  variance_map.end() &&
	 frb_event.event_ts.find(*ts) !=  frb_event.event_ts.end() &&
	 var > var_low_cutoff && 
	 var < var_high_cutoff){

      /**
      G3Timestream ts_local(frb_event.event_ts.at(*ts));
      //mean filter and take absolute value
      double mean_val = 0;
      for (size_t i=0; i < ts_local.size(); i++) mean_val += ts_local[i];
      mean_val /= (double) ts_local.size();
      for (size_t i=0; i < ts_local.size(); i++) {
	ts_local[i] = fabs( ts_local[i]);
      }
      sum_ts = sum_ts + ts_local;
      **/

      sum_ts = sum_ts + *frb_event.event_ts.at(*ts);
      sum_variance += var;
    }
  }
  if (sum_variance == 0){
	  return 0;
  }
  double model_ll =  get_individual_delta_heavi_poly_fit_ll( sum_ts, sum_variance, ev_index,
							     poly_order, fit_model_len,
							     true, include_heavi);

  double base_ll =  get_individual_delta_heavi_poly_fit_ll( sum_ts, sum_variance, ev_index,
							    poly_order, fit_model_len,
							    false, include_heavi);
  return model_ll - base_ll;
}


G3MapVectorStringPtr get_groups_of_importance( const G3MapString & bolo_to_squid,
					       const G3MapString & bolo_to_wafer,
					       const G3MapVectorString & squid_to_bolo,
					       const G3MapVectorString & wafer_to_bolo,
					       const FrbEventInfo & ev
					       ){
  log_debug("entering");

  std::unordered_set<std::string> unique_squids;
  std::unordered_set<std::string> unique_wafers;
  for (auto it = ev.triggers.begin(); it != ev.triggers.end(); it++){
    unique_squids.insert(bolo_to_squid.at(*it));
    unique_wafers.insert(bolo_to_wafer.at(*it));
  }
  G3MapVectorStringPtr groups(new G3MapVectorString);
  for (auto it = unique_squids.begin(); it != unique_squids.end(); it++){
    (*groups)[*it] = squid_to_bolo.at(*it);
  }
  for (auto it = unique_wafers.begin(); it != unique_wafers.end(); it++){
    (*groups)[*it] = wafer_to_bolo.at(*it);
  }
  log_debug("exiting");

  return groups;
}





void GroupLoglikelihoodFilter::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("entering");
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_in_key_);
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    G3MapStringConstPtr bolo_to_squid = frame->Get<G3MapString>(bolo_to_squid_key_);
    G3MapStringConstPtr bolo_to_wafer = frame->Get<G3MapString>(bolo_to_wafer_key_);
    G3MapVectorStringConstPtr squid_to_bolo = frame->Get<G3MapVectorString>(squid_to_bolo_key_);
    G3MapVectorStringConstPtr wafer_to_bolo = frame->Get<G3MapVectorString>(wafer_to_bolo_key_);
    G3MapDoubleConstPtr variance_map = frame->Get<G3MapDouble>(variance_map_key_);

    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      FrbEventInfo new_ev = *ev;//needs to be copied since debug info it being added
      G3MapVectorStringPtr groups_of_import = get_groups_of_importance(*bolo_to_squid,
								       *bolo_to_wafer,
								       *squid_to_bolo,
								       *wafer_to_bolo,
								       *ev );


      bool boot_the_ev = false;
      for (auto grp = groups_of_import->begin(); grp != groups_of_import->end(); grp++){
	double ll_diff =  get_group_sum_loglikelihood(new_ev, grp->second,
						      *variance_map,
						      true,
						      11, 1, true,
						      variance_low_cutoff_,
						      variance_high_cutoff_ );
	new_ev.debug_info[grp->first] = ll_diff;

	if (ll_diff > loglikelihood_cutoff_ || ll_diff == 0){
	  boot_the_ev = true;
	}
      }
      if (! boot_the_ev){
	out_lst->push_back(new_ev);
      }
    }
    frame->Put(frb_out_key_, out_lst);
  }
  out.push_back(frame);
}






void FilterEventsWithOtherModelSize::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_ );
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      bool filter_ev = false;
      for (auto it = ev->triggers.begin(); it != ev->triggers.end(); it++){
	const G3Timestream & ts = *ev->event_ts.at(*it);
	size_t ev_len = ts.size();
	//g3_assert(ev_len%2 == 1); // i don't know if this will line up right with even lengths so might as well fail than lie
	size_t event_ind = ev_len/2;
	double variance = frame->Get<G3MapDouble>(variance_key_)->at(*it);
	double mll = get_individual_delta_heavi_poly_fit_ll(ts, variance,
							    event_ind,
							    1, new_model_size_,
							    true, true); //delta, heavi

	double bll = get_individual_delta_heavi_poly_fit_ll(ts, variance,
							    event_ind,
							    1, new_model_size_,
							    false, true); //delta, heavi

	double delta = mll-bll;
	if (delta < loglike_cutoff_){
	  filter_ev = true;
	  break;
	}
      }
      if (!filter_ev) out_lst->push_back(*ev);
    }
    frame->Put( out_frb_event_key_, out_lst);
  }
  out.push_back(frame);
}






SmartAmpFilter::SmartAmpFilter(std::string frb_event_key,
			       std::string frb_event_out_key,
			       std::string variance_path, 
			       double n_sigma_away,
			       double max_amp_err,
			       bool invert_result
			       ) :
  frb_event_key_(frb_event_key), frb_event_out_key_(frb_event_out_key),
  variance_path_(variance_path), n_sigma_away_(n_sigma_away),max_amp_err_(max_amp_err),
  invert_result_(invert_result)
{}


void SmartAmpFilter::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("Enter");
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_ );
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    G3MapDoubleConstPtr variances = frame->Get<G3MapDouble>(variance_path_);
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      g3_assert(ev->triggers.size() > 0);
      g3_assert(ev->is_valid());
      size_t big_ind = 0;
      double big_size = ev->trigger_ll_model.at(0) - ev->trigger_ll_baseline.at(0);

      for (size_t j=1; j < ev->triggers.size(); j++){
	double ll =ev->trigger_ll_model.at(j) - ev->trigger_ll_baseline.at(j);
	if (ll > big_size){
	  big_size = ll;
	  big_ind = j;
	}
      }
      bool include_ev = true;

      log_debug("variance grab");
      double bigvar = variances->at(ev->triggers.at(big_ind));
      log_debug("amp grab");
      double bigamp = ev->trigger_amplitude.at(big_ind);
      log_debug("loop");

      for (size_t j=0; j < ev->triggers.size(); j++){
	if (j == big_ind) continue;
	double var = variances->at(ev->triggers.at(j));
	double amp =  ev->trigger_amplitude.at(j);
	
	double var_delta = sqrt(bigvar + var) * n_sigma_away_;
	double amp_delta = bigamp * max_amp_err_;
	double delta = MAX(var_delta, amp_delta);

	if (fabs(bigamp-amp) > delta){
	  include_ev = false;
	  break;
	}
      }
      if (invert_result_ ^ include_ev) { 
	      out_lst->push_back(*ev);
      } 
    }
    frame->Put( frb_event_out_key_, out_lst);
  }
  out.push_back(frame);

}




CutYDetectorEvents::CutYDetectorEvents(std::string frb_event_key,
				       std::string out_frb_event_key,
				       bool cut_x_instead
				       ) : 
  frb_event_key_(frb_event_key),  out_frb_event_key_(out_frb_event_key),
  cut_x_instead_(cut_x_instead)
{}

void CutYDetectorEvents::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("Entering process");
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_);
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      log_debug("ev");
      g3_assert(ev->triggers.size() == ev->trigger_ll_model.size());
      g3_assert(ev->triggers.size() == ev->trigger_ll_baseline.size());
      g3_assert(ev->triggers.size() == ev->trigger_amplitude.size());

      char match_char = cut_x_instead_ ? 'Y' : 'X';
      
      FrbEventInfo info;
      info.copy_scalars(*ev);
      for (size_t i=0; i < ev->triggers.size(); i++){
	std::string det_name = ev->triggers[i];
	if ( det_name[det_name.size()-1] == match_char ){
	  info.triggers.push_back(ev->triggers[i]);
	  info.trigger_ll_model.push_back(ev->trigger_ll_model[i]);
	  info.trigger_ll_baseline.push_back(ev->trigger_ll_baseline[i]);
	  info.trigger_amplitude.push_back(ev->trigger_amplitude[i]);
	  info.trigger_heavi_amp.push_back(ev->trigger_heavi_amp[i]);
	}
      }
      if (info.triggers.size() > 0) out_lst->push_back(info);      
    }
    frame->Put(out_frb_event_key_, out_lst);
  }
  out.push_back(frame);
}



FilterYDetectorEvents::FilterYDetectorEvents(std::string frb_event_key,
				       std::string out_frb_event_key,
				       bool cut_x_instead
				       ) : 
  frb_event_key_(frb_event_key),  out_frb_event_key_(out_frb_event_key),
  cut_x_instead_(cut_x_instead)
{}

void FilterYDetectorEvents::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("Entering process");
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_);
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      log_debug("ev");
      g3_assert(ev->triggers.size() == ev->trigger_ll_model.size());
      g3_assert(ev->triggers.size() == ev->trigger_ll_baseline.size());
      g3_assert(ev->triggers.size() == ev->trigger_amplitude.size());

      char match_char = cut_x_instead_ ? 'Y' : 'X';
      bool include = true;
      for (size_t i=0; i < ev->triggers.size(); i++){
	std::string det_name = ev->triggers[i];
	if ( det_name[det_name.size()-1] != match_char ){
	  include = false;
	}
      }
      if (include) out_lst->push_back(*ev);      
    }
    frame->Put(out_frb_event_key_, out_lst);
  }
  out.push_back(frame);
}




CutGlitchyDetectorsFrb::CutGlitchyDetectorsFrb(std::string frb_event_key, 
					       std::string observation_name_key,
					       double prob_cutoff,
					       int minimum_for_cut,
					       std::string new_frb_event_key) :
  frb_event_key_(frb_event_key), 
  new_frb_event_key_(new_frb_event_key),
  observation_name_key_(observation_name_key),
  prob_cutoff_(prob_cutoff), minimum_for_cut_(minimum_for_cut){
  observation_name_ = "";
  events_ = G3VectorFrbEventInfoPtr(new G3VectorFrbEventInfo);
  
}

void CutGlitchyDetectorsFrb::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("Enter");
  if (frame->type == G3Frame::Scan && observation_name_ == ""){
    observation_name_ = frame->Get<G3String>(observation_name_key_)->value;
  }

  if (frame->type == G3Frame::EndProcessing || 
      (frame->type == G3Frame::Scan && observation_name_ != frame->Get<G3String>(observation_name_key_)->value )
      ){

    if (frame->type == G3Frame::Scan) observation_name_ = frame->Get<G3String>(observation_name_key_)->value;

    std::unordered_set<std::string> bad_det;
    get_too_many_frb_events(events_,  prob_cutoff_, minimum_for_cut_, bad_det);
    for (auto iframe = collected_frames_.begin(); iframe != collected_frames_.end(); iframe++ ) {
      G3VectorFrbEventInfoConstPtr evs = (**iframe).Get<G3VectorFrbEventInfo>(frb_event_key_);
      G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());

      for (auto ev = evs->begin(); ev != evs->end(); ev++){
	FrbEventInfo info;
	info.copy_scalars(*ev);
	for (size_t i=0; i < ev->triggers.size(); i++){
	  if (bad_det.count(ev->triggers[i]) == 0 ){
	    info.triggers.push_back(ev->triggers[i]);
	    info.trigger_ll_model.push_back(ev->trigger_ll_model[i]);
	    info.trigger_ll_baseline.push_back(ev->trigger_ll_baseline[i]);
	    info.trigger_amplitude.push_back(ev->trigger_amplitude[i]);
	    info.trigger_heavi_amp.push_back(ev->trigger_heavi_amp[i]);
	  }
	}
	if (info.triggers.size() > 0) out_lst->push_back(info);
      }
      (*iframe)->Put(new_frb_event_key_, out_lst);
      out.push_back( *iframe);
    }
    collected_frames_.clear();
    events_->clear();
  } 

  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_);
    (*events_).insert((*events_).end(), (*evs).begin(), (*evs).end());
    collected_frames_.push_back(frame);
  }else{
    out.push_back(frame);
  }
}

CutNegAmpInEvents::CutNegAmpInEvents(std::string frb_event_key,
				     std::string out_frb_event_key ) : 
  frb_event_key_(frb_event_key),  out_frb_event_key_(out_frb_event_key)
{}

void CutNegAmpInEvents::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("Entering process");
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_);
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      log_debug("ev");
      g3_assert(ev->triggers.size() == ev->trigger_ll_model.size());
      g3_assert(ev->triggers.size() == ev->trigger_ll_baseline.size());
      g3_assert(ev->triggers.size() == ev->trigger_amplitude.size());
      
      int n_neg=0;
      int n_pos=0;
      for (size_t i=0; i < ev->triggers.size(); i++){
	if ( ev->trigger_amplitude[i] <= 0  ){
	  n_neg++;
	} else{
	  n_pos++;
	}
      }
      if (n_neg == 0){
	out_lst->push_back(*ev);
	continue;
      }
      FrbEventInfo info;
      info.copy_scalars(*ev);
      info.triggers = std::vector<std::string>(n_pos);
      info.trigger_ll_model = std::vector<double>(n_pos);
      info.trigger_ll_baseline = std::vector<double>(n_pos);
      info.trigger_amplitude = std::vector<double>(n_pos);
      info.trigger_heavi_amp = std::vector<double>(n_pos);
      int j_tmp = 0;
      for (size_t i=0; i < ev->triggers.size(); i++){
	if ( ev->trigger_amplitude[i] > 0  ){
	  info.triggers[j_tmp] = ev->triggers[i];
	  info.trigger_ll_model[j_tmp] = ev->trigger_ll_model[i];
	  info.trigger_ll_baseline[j_tmp] = ev->trigger_ll_baseline[i];
	  info.trigger_amplitude[j_tmp] = ev->trigger_amplitude[i];
	  info.trigger_heavi_amp[j_tmp] = ev->trigger_heavi_amp[i];
	  j_tmp++;
	}
      }
      if (info.triggers.size() > 0){
	out_lst->push_back(info);      
      }
    }
    frame->Put(out_frb_event_key_, out_lst);
  }
  out.push_back(frame);
}





CutAmpInEvents::CutAmpInEvents(std::string frb_event_key,
			       std::string out_frb_event_key,
			       double min_amp, double max_amp
			       ) : 
  frb_event_key_(frb_event_key),  out_frb_event_key_(out_frb_event_key),
  min_amp_(min_amp), max_amp_(max_amp)
{}

void CutAmpInEvents::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("Entering process");
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_);
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      log_debug("ev");
      g3_assert(ev->triggers.size() == ev->trigger_ll_model.size());
      g3_assert(ev->triggers.size() == ev->trigger_ll_baseline.size());
      g3_assert(ev->triggers.size() == ev->trigger_amplitude.size());

      FrbEventInfo info;
      info.copy_scalars(*ev);
      for (size_t i=0; i < ev->triggers.size(); i++){
	//log_warn("amp %lf\n", ev->trigger_amplitude[i]);
	if ( ev->trigger_amplitude[i] > min_amp_ && ev->trigger_amplitude[i] < max_amp_  ){
	  info.triggers.push_back(ev->triggers[i]);
	  info.trigger_ll_model.push_back(ev->trigger_ll_model[i]);
	  info.trigger_ll_baseline.push_back(ev->trigger_ll_baseline[i]);
	  info.trigger_amplitude.push_back(ev->trigger_amplitude[i]);
	  info.trigger_heavi_amp.push_back(ev->trigger_heavi_amp[i]);
	}
      }
      if (info.triggers.size() > 0) out_lst->push_back(info);      
    }
    frame->Put(out_frb_event_key_, out_lst);
  }
  out.push_back(frame);
}



FilterForDetectorPair::FilterForDetectorPair(std::string frb_event_key,
					     std::string frb_event_out_key) :
  frb_event_key_(frb_event_key),
  frb_event_out_key_(frb_event_out_key)
{}


void FilterForDetectorPair::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("Enter");
  if (frame->type == G3Frame::Scan){
    g3_assert(frame->Has(frb_event_key_));
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_ );
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      if (ev->triggers.size() == 2 && 
	  ( ev->triggers[0].substr(0, ev->triggers[0].size()-1) == 
	    ev->triggers[1].substr(0, ev->triggers[1].size()-1))){

	log_notice("Found Pair %s %s", ev->triggers[0].c_str(), ev->triggers[1].c_str());
	out_lst->push_back(*ev);
      }
    }
    frame->Put( frb_event_out_key_, out_lst);
  }
  out.push_back(frame);
}



void FilterEventsWithBigNegDeviation::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_ );
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      bool filter_ev = false;
      for (auto it = ev->triggers.begin(); it != ev->triggers.end(); it++){
	G3Timestream ts = mean_filter_timestream(*ev->event_ts.at(*it));
	if (fabs( get_min(ts)) > fabs( get_max(ts))){
	  filter_ev = true;
	  break;
	}
      }
      if (!filter_ev) out_lst->push_back(*ev);
    }
    frame->Put( out_frb_event_key_, out_lst);
  }
  out.push_back(frame);
}





FilterNumDets::FilterNumDets(std::string frb_event_key,
			     std::string frb_event_out_key,
			     size_t n_dets) :
  frb_event_key_(frb_event_key), frb_event_out_key_(frb_event_out_key),
  n_dets_(n_dets)
{}


void FilterNumDets::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
  log_debug("Enter");
  if (frame->type == G3Frame::Scan){
    G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_ );
    G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
    for (auto ev = evs->begin(); ev != evs->end(); ev++){
      if (ev->triggers.size() == n_dets_){
	out_lst->push_back(*ev);
      }
    }
    frame->Put( frb_event_out_key_, out_lst);
  }
  out.push_back(frame);
}












void MakeOtherThresholds::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	log_debug("Entering process");
	if (frame->type == G3Frame::Scan){
		G3VectorFrbEventInfoConstPtr evs = frame->Get<G3VectorFrbEventInfo>(frb_event_key_);
		G3VectorFrbEventInfoPtr out_lst(new G3VectorFrbEventInfo());
		for (auto ev = evs->begin(); ev != evs->end(); ev++){
			log_debug("ev");
			g3_assert(ev->triggers.size() == ev->trigger_ll_model.size());
			g3_assert(ev->triggers.size() == ev->trigger_ll_baseline.size());
			g3_assert(ev->triggers.size() == ev->trigger_amplitude.size());

			FrbEventInfo info;
			info.copy_scalars(*ev);
			/**
			info.triggers = std::vector<std::string>();
			info.trigger_ll_model = std::vector<double>();
			info.trigger_ll_baseline = std::vector<double>();
			info.trigger_amplitude = std::vector<double>();

			int j_tmp = 0;
			**/
			double max_event_ll = 0;
			for (size_t i=0; i < ev->triggers.size(); i++){
				double delta = (ev->trigger_ll_model[i] - ev->trigger_ll_baseline[i]);
 				if ( delta > other_threshold_ ){
					max_event_ll = delta > max_event_ll ? delta : max_event_ll;
					info.triggers.push_back(ev->triggers[i]);
					info.trigger_ll_model.push_back(ev->trigger_ll_model[i]);
					info.trigger_ll_baseline.push_back(ev->trigger_ll_baseline[i]);
					info.trigger_amplitude.push_back(ev->trigger_amplitude[i]);
					info.trigger_heavi_amp.push_back(ev->trigger_heavi_amp[i]);
				}
			}
			if (max_event_ll > event_threshold_) {
				out_lst->push_back(info);      
			}
		}
		frame->Put(frb_event_out_key_, out_lst);
	}
	out.push_back(frame);
	
}
