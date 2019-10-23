#pragma once

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <boost/math/distributions/poisson.hpp>
#include <memory>
#include <list>
#include <deque>
#include <unordered_set>

#include <G3Timestream.h>
#include <G3Module.h>
#include <G3Vector.h>
#include <G3Map.h>
#include <G3Logging.h>
#include <G3TimeStamp.h>

#include <todfilter/fftutils.h>
#include <todfilter/polyutils.h>

#include <core/pybindings.h>
#include <core/serialization.h>

namespace bp = boost::python;

struct FrbDetInfo : public G3FrameObject{
	std::string bid;
	double significance;
	
	double ll_model;
	double ll_baseline;
	
	double amplitude;
	double heavi_amp;

	double q_significance;
	double q_ll_model;
	double q_ll_baseline;
	double q_amplitude;
	double q_heavi_amp;
	
	int board_id;
	int module_id;	
	double band;

	double squid_sig;
	double squid_sig_max;
	double wafer_sig;
	
	double ra;
	double dec;
	
	bool was_injected;

	size_t scan_index;
	
	bool operator==(const FrbDetInfo & other){
		return (bid == other.bid &&
			significance == other.significance &&
			ll_model == other.ll_model &&
			ll_baseline == other.ll_baseline &&
			amplitude == other.amplitude &&
			heavi_amp == other.heavi_amp &&
			board_id == other.board_id &&
			module_id == other.module_id &&
			squid_sig == other.squid_sig &&
			wafer_sig == other.wafer_sig &&
			ra == other.ra &&
			dec == other.dec
			);
	}
	template <class A> void serialize(A &ar, unsigned v) {
		using namespace cereal;

		G3_CHECK_VERSION(v);
		ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
		ar & make_nvp("bid", bid);
		
		ar & make_nvp("significance", significance);

		ar & make_nvp("ll_model", ll_model);
		ar & make_nvp("ll_baseline", ll_baseline);

		ar & make_nvp("amplitude", amplitude);
		ar & make_nvp("heavi_amp", heavi_amp);

		ar & make_nvp("board_id", board_id);
		ar & make_nvp("module_id", module_id);
		ar & make_nvp("band", band);

		ar & make_nvp("squid_sig", squid_sig);
		ar & make_nvp("wafer_sig", wafer_sig);

		ar & make_nvp("ra", ra);
		ar & make_nvp("dec", dec);

		ar & make_nvp("was_injected", was_injected);
		
		if ( v > 4) {
			ar & make_nvp("scan_index", scan_index);
		}
		if (v > 5) {
			ar & make_nvp("q_significance", q_significance);
			ar & make_nvp("q_ll_model", q_ll_model);
			ar & make_nvp("q_ll_baseline", q_ll_baseline);
			ar & make_nvp("q_amplitude", q_amplitude);
			ar & make_nvp("q_heavi_amp", q_heavi_amp);
		}
		if (v>6) {
			ar & make_nvp("squid_sig_max", squid_sig_max);			
		}
	}
};

G3_POINTERS(FrbDetInfo);
G3_SERIALIZABLE(FrbDetInfo,7);
G3VECTOR_OF(FrbDetInfo, G3VectorFrbDetInfo);

struct FrbEventInfo : public G3FrameObject{
	std::vector<FrbDetInfo> det_info;

	std::string observation_name;
	int64_t observation_number;
	size_t scan_number;
	size_t scan_index;

	G3Time event_time;

	void filter_det_evs( std::vector<int> is_bad) {
		if (is_bad.size() == 0) return;
		g3_assert(is_bad.size() == det_info.size());
		size_t n_bad = 0;
		size_t tsize = det_info.size();
		for (int i = is_bad.size() - 1; i >= 0; i--) {
			if (is_bad[i]) {
				n_bad++;
				g3_assert(tsize - n_bad >= 0);
				det_info[i] = det_info[tsize - n_bad];
			}
		}
		det_info.resize(det_info.size() - n_bad);
	}

	template <class A> void serialize(A &ar, unsigned v) {
		using namespace cereal;
		G3_CHECK_VERSION(v);
		ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));

		ar & make_nvp("det_info", det_info);
		ar & make_nvp("observation_name",observation_name);
		ar & make_nvp("observation_number",observation_number);
		ar & make_nvp("scan_number",scan_number);
		ar & make_nvp("scan_index",scan_index);


		if (v > 6) {
			ar & make_nvp("event_time", event_time);
		}
	}
	
	bool operator==(const FrbEventInfo & other){
		return (other.observation_name == observation_name &&
			other.observation_number == observation_number &&
			other.scan_number == scan_number &&
			other.scan_index == scan_index);
	}
};

G3_POINTERS(FrbEventInfo);
G3_SERIALIZABLE(FrbEventInfo,7);
G3VECTOR_OF(FrbEventInfo, G3VectorFrbEventInfo);



class PolyLikelihoodFiller : public G3Module{
public:
	PolyLikelihoodFiller(size_t model_len, size_t poly_order, 
			     bool include_heaviside, 
			     bool include_delta, 
			     std::string ts_key, 
			     std::string variance_key,
			     std::string loglike_output_key,
			     std::string amp_map_output_key,
			     std::string hamp_map_output_key
		);
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	size_t model_len_;
	size_t poly_order_;
	bool include_heaviside_;
	bool include_delta_;
	std::string ts_key_;
	std::string variance_key_;
	std::string loglike_output_key_;
	std::string amp_map_output_key_;
	std::string hamp_map_output_key_;
};

class DeltaEventHunter : public G3Module{
public:
	DeltaEventHunter(std::string ll_model_key,  
			 std::string ll_base_key,  
			 double trigger_thresh, 
			 double other_det_thresh,
			 size_t min_distance,
			 std::string output_event_key,
			 std::string fit_amp_key, 
			 std::string fit_hamp_key, 
			 int search_width
		   );
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	std::string ll_model_key_;  
	std::string ll_base_key_;  
	double trigger_thresh_; 
	double other_det_thresh_;
	size_t min_distance_;
	std::string output_event_key_;
	std::string fit_amp_key_;
	std::string fit_hamp_key_;
	int search_width_;
};





