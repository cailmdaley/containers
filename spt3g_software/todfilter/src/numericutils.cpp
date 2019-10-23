#include <pybindings.h>

#include <todfilter/numericutils.h>
#include <G3Frame.h>
#include <G3Timestream.h>
#include <G3Map.h>
#include <G3Vector.h>

#include <vector>
#include <math.h>

#ifdef OPENMP_FOUND
#include <omp.h>
#endif



void average_psd_bins(double freq_sep,
		  const G3MapVectorDouble & psds,
		  const G3VectorDouble & lower_bounds,
		  const G3VectorDouble & upper_bounds,
		  G3MapVectorDouble & grouped_sum_vals
	) {
	g3_assert(lower_bounds.size() == upper_bounds.size());

	grouped_sum_vals.clear();
	for (auto it=psds.begin(); it != psds.end(); it++){
		const std::vector<double> & psd = it->second;
		G3VectorDouble sum_vals(lower_bounds.size(), 0);
		for (size_t bounds_i = 0; bounds_i < lower_bounds.size(); bounds_i ++){
			double lower_bound = lower_bounds[bounds_i];
			double upper_bound = upper_bounds[bounds_i];
			size_t lower_ind = (size_t) ceilf( lower_bound / freq_sep);
			size_t upper_ind = (size_t) floorf( upper_bound / freq_sep);
			
			if (upper_ind < lower_ind) continue;
			double n_bins = upper_ind - lower_ind + 1;
			for (size_t j=lower_ind; j <= upper_ind; j++){
				sum_vals[bounds_i] += psd[j];
			}
			
			sum_vals[bounds_i] /= (n_bins);
		}
		grouped_sum_vals[it->first] = sum_vals;
	}
}

bool is_power_of_two(unsigned int x) {
	return ((x != 0) && ((x & (~x + 1)) == x));
}

void sum_reduce_vec_to_size(double * iv, size_t iv_size, size_t new_size){
	//deprecated?
	if (iv_size == new_size) return;
	if (iv_size % new_size != 0){
		log_fatal("sum reduce can only handle things that evenly divide in %zu %zu", 
			  iv_size, new_size);
	}
	size_t nper = iv_size / new_size;

	G3VectorDouble ov_tmp(new_size, 0);
	for (int i=0; i < new_size; i++) {
		for (int j=0; j < nper; j++) {
			ov_tmp[i] += iv[i * nper + j];
		}
	}
	for (int i=0; i < new_size; i++) 
		iv[i] = ov_tmp[i];
}


double get_lower_median(G3Timestream v) {
	return quick_select<double>(&(v[0]), v.size());
}


double get_mad_to_variance_correction_normal_dist(){
	return 1.4826;
}


G3Timestream mean_filter_timestream(const G3Timestream & ts){
	G3Timestream out_ts(ts);
	double mean_val = 0;
	for (size_t i=0; i < out_ts.size(); i++) mean_val += out_ts[i];
	mean_val /= (double) ts.size();
	for (size_t i=0; i < out_ts.size(); i++) out_ts[i] /= mean_val;
	return out_ts;
}

double get_mean(const G3Timestream & ts){
	double mean_val = 0;
	for (size_t i=0; i < ts.size(); i++) mean_val += ts.at(i);
	mean_val /= (double) ts.size();
	return mean_val;
}

double get_min(const G3Timestream & ts){
	double min_val = ts[0];
	for (size_t i=1; i < ts.size(); i++) min_val =  ts[i] < min_val ? ts[i] : min_val;
	return min_val;
}

double get_max(const G3Timestream & ts){
	double max_val = ts[0];
	for (size_t i=1; i < ts.size(); i++) max_val =  ts[i] > max_val ? ts[i] : max_val;
	return max_val;
}

double get_mad_std(G3TimestreamConstPtr ts) {
	return get_mad_to_variance_correction_normal_dist() * 
		compute_median_absolute_deviation<double>(&((*ts)[0]), ts->size());
}

template void find_local_maxima<double>(const std::vector<double> & in_vec, 
					int min_index, size_t max_index,
					int & out_index, double & out_val);



MadVarianceAdder::MadVarianceAdder(std::string ts_key, std::string variance_output_key) :
	variance_output_key_(variance_output_key), ts_key_(ts_key){}

void MadVarianceAdder::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	if (frame->type == G3Frame::Scan){
		const double mad_corr = get_mad_to_variance_correction_normal_dist();
		G3TimestreamMapConstPtr ts_map = frame->Get<G3TimestreamMap>(ts_key_);
		G3MapDoublePtr mad_map(new G3MapDouble);
		for (auto it=ts_map->begin(); it != ts_map->end(); it++){
			double mad = mad_corr * compute_median_absolute_deviation<double>(&((*it->second)[0] ), it->second->size());
			(*mad_map)[it->first] = mad * mad;
		}
		frame->Put(variance_output_key_, mad_map);
	}
	out.push_back(frame);
}


VarianceAdder::VarianceAdder(std::string ts_key, std::string variance_output_key) :
	variance_output_key_(variance_output_key), ts_key_(ts_key){}

void VarianceAdder::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	if (frame->type == G3Frame::Scan){
		G3TimestreamMapConstPtr ts_map = frame->Get<G3TimestreamMap>(ts_key_);
		G3MapDoublePtr var_map(new G3MapDouble);
		for (auto it=ts_map->begin(); it != ts_map->end(); it++){
			double var =0;
			double mean = 0;
			for (size_t i=0; i < it->second->size(); i++){
				mean += (*it->second)[i];
				var += (*it->second)[i] * (*it->second)[i];
			}
			mean /= it->second->size();
			var /= it->second->size();
			var = var - mean*mean;
			(*var_map)[it->first] = var;
		}
		frame->Put(variance_output_key_, var_map);
	}
	out.push_back(frame);
}

namespace bp = boost::python;

PYBINDINGS("todfilter") {
	bp::class_<MadVarianceAdder, bp::bases<G3Module>,
	           boost::shared_ptr<MadVarianceAdder>,
	           boost::noncopyable>("MadVarianceAdder",
	               bp::init<std::string, std::string>
	               (MAD_VARIANCE_ADDER_DOCSTR,
	                   (bp::arg("ts_key"), bp::arg("variance_output_key"))
	               ))
	.def("Process", &MadVarianceAdder::Process);

	bp::implicitly_convertible< boost::shared_ptr<MadVarianceAdder>, G3ModulePtr>();

	bp::class_<VarianceAdder, bp::bases<G3Module>,
	           boost::shared_ptr<VarianceAdder>,
	           boost::noncopyable>("VarianceAdder",
	               bp::init<std::string, std::string>
	               (VARIANCE_ADDER_DOCSTR,
	                   (bp::arg("ts_key"), bp::arg("variance_output_key"))
	               ))
	.def("Process", &VarianceAdder::Process);

	bp::implicitly_convertible< boost::shared_ptr<VarianceAdder>, G3ModulePtr>();

	bp::def("get_mad_std", get_mad_std, 
		bp::arg("timestream"), 
		"Estimates the standard deviation of a timestream with a median absolute "
		"deviation estimator.  This is an outlier insensitive estimate that "
		"makes the assumption that the noise is Gaussian.");
	bp::def("get_lower_median", 
		get_lower_median, bp::arg("timestream"), 
		"Returns the median of the data.  If there are an even number of samples, "
		"it doesn't average the two middle values and just returns the smaller of "
		"the two value.");
	bp::def("average_psd_bins", average_psd_bins,
		(bp::arg("freq_sep"), bp::arg("psds"), bp::arg("lower_bounds"),
		 bp::arg("upper_bounds"), bp::arg("grouped_sum_vals")),
		"Averages the psd values over the bins specified by lower_bounds and upper_bounds.  Returns teh value in grouped_sum_vals."
		);

	bp::def("get_mean", get_mean);
}
