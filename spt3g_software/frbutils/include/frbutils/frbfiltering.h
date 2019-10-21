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

#include <todfilter/polyutils.h>
#include <todfilter/fftutils.h>
#include "frbutils.h"

#include <core/pybindings.h>
namespace bp = boost::python;


class GroupLoglikelihoodFilter : public G3Module { 
 public:
  GroupLoglikelihoodFilter(std::string frb_in_key,
			   std::string frb_out_key,
			   std::string bolo_to_squid_key,
			   std::string bolo_to_wafer_key,
			   std::string squid_to_bolo_key,
			   std::string wafer_to_bolo_key,
			   std::string variance_map_key,
			   double loglikelihood_cutoff,
			   double variance_low_cutoff,
			   double variance_high_cutoff
			   ) :
  frb_in_key_(frb_in_key), frb_out_key_(frb_out_key),
    bolo_to_squid_key_(bolo_to_squid_key),
    bolo_to_wafer_key_(bolo_to_wafer_key),
    squid_to_bolo_key_(squid_to_bolo_key),
    wafer_to_bolo_key_(wafer_to_bolo_key),
    variance_map_key_(variance_map_key),
    loglikelihood_cutoff_(loglikelihood_cutoff),
    variance_low_cutoff_(variance_low_cutoff),
    variance_high_cutoff_(variance_high_cutoff)
    
  {}
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_in_key_;
  std::string frb_out_key_;

  std::string bolo_to_squid_key_;
  std::string bolo_to_wafer_key_;
  std::string squid_to_bolo_key_;
  std::string wafer_to_bolo_key_;
  std::string variance_map_key_;
  double loglikelihood_cutoff_;
  double variance_low_cutoff_;
  double variance_high_cutoff_;
};


EXPORT_G3MODULE("frbutils", 
		GroupLoglikelihoodFilter, 
		(bp::init<std::string, std::string, std::string, 
		   std::string, std::string, std::string, std::string, double,
		   double, double >(
		(bp::arg("frb_in_key"),
				      bp::arg("frb_out_key"),
				      bp::arg("bolo_to_squid_key"),
				      bp::arg("bolo_to_wafer_key"),
				      bp::arg("squid_to_bolo_key"),
				      bp::arg("wafer_to_bolo_key"),
				      bp::arg("variance_map_key"),
				      bp::arg("loglikelihood_cutoff"),
				      bp::arg("variance_low_cutoff"),
				      bp::arg("variance_high_cutoff")
		 ))),
		"Hands are weird");



class FilterEventsWithOtherModelSize : public G3Module { 
 public:
 FilterEventsWithOtherModelSize(std::string frb_event_key, 
				std::string out_frb_event_key,
				std::string variance_key,
				size_t new_model_size,
				double loglike_cutoff) :
  frb_event_key_(frb_event_key), out_frb_event_key_(out_frb_event_key),
    new_model_size_(new_model_size), loglike_cutoff_(loglike_cutoff),
    variance_key_(variance_key)
  {}

  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string out_frb_event_key_;
  size_t new_model_size_;
  double loglike_cutoff_;
  std::string variance_key_;
};

EXPORT_G3MODULE("frbutils", 
		FilterEventsWithOtherModelSize, 
		(init<std::string, std::string, std::string, size_t, double>(
		(bp::arg("frb_event_key"),
		 bp::arg("out_frb_event_key"),
		 bp::arg("variance_key"),
		 bp::arg("new_model_size"),
		 bp::arg("loglike_cutoff")
		 ))),
		"Hands are weird");



class CutGlitchyDetectorsFrb : public G3Module{
 public:
  CutGlitchyDetectorsFrb( std::string frb_event_key, 
			  std::string observation_name_key, 
			  double prob_cutoff,
			  int minimum_for_cut,
			  std::string new_frb_event_key);
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string new_frb_event_key_;
  std::string observation_name_key_;  
  std::string observation_name_;
  double prob_cutoff_;
  int minimum_for_cut_;
  G3VectorFrbEventInfoPtr events_;
  std::deque<G3FramePtr> collected_frames_;
};


EXPORT_G3MODULE("frbutils", 
		CutGlitchyDetectorsFrb, 
		( bp::init<std::string, std::string,
		  double, int, std::string>(
		(bp::arg("frb_event_key"),
		 bp::arg("observation_name_key"),
		 bp::arg("prob_cutoff"),
		 bp::arg("minimum_for_cut"),
		 bp::arg("new_frb_event_key")
		 ))),
		"Hands are weird");




class CutNegAmpInEvents : public G3Module { 
 public:
  CutNegAmpInEvents(std::string frb_event_key, 
		    std::string out_frb_event_key);
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string out_frb_event_key_;
};


EXPORT_G3MODULE("frbutils", 
		CutNegAmpInEvents,
		(
		 bp::init<std::string, std::string>
		 ( ( bp::arg("frb_event_key"),
		     bp::arg("out_frb_event_key")
		     ))
		 ),
		"Hands are weird");



class CutAmpInEvents : public G3Module { 
 public:
  CutAmpInEvents(std::string frb_event_key, 
		    std::string out_frb_event_key,
		    double min_amp, double max_amp
		    );
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string out_frb_event_key_;
  double min_amp_;
  double max_amp_;
};

EXPORT_G3MODULE("frbutils", 
		CutAmpInEvents,
		(
		 bp::init<std::string, std::string, double, double>
		 ( ( bp::arg("frb_event_key"),
		     bp::arg("out_frb_event_key"),
		     bp::arg("min_amp"),
		     bp::arg("max_amp")
		     )
		   )
		 ),
		"Hands are weird");




class CutYDetectorEvents : public G3Module { 
 public:
  CutYDetectorEvents(std::string frb_event_key, 
		     std::string out_frb_event_key,
		     bool cut_x_instead
		     );
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string out_frb_event_key_;
  bool cut_x_instead_;
};


EXPORT_G3MODULE("frbutils", 
		CutYDetectorEvents,
		(
		 bp::init<std::string, std::string, bool>
		 ( ( bp::arg("frb_event_key"),
		     bp::arg("out_frb_event_key"),
		     bp::arg("cut_x_instead")
		     )
		   )
		 ),
		"Hands are weird");




class FilterYDetectorEvents : public G3Module { 
 public:
  FilterYDetectorEvents(std::string frb_event_key, 
		     std::string out_frb_event_key,
		     bool cut_x_instead
		     );
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string out_frb_event_key_;
  bool cut_x_instead_;
};

EXPORT_G3MODULE("frbutils", 
		FilterYDetectorEvents,
		(
		 bp::init<std::string, std::string, bool>
		 ( ( bp::arg("frb_event_key"),
		     bp::arg("out_frb_event_key"),
		     bp::arg("cut_x_instead")
		     )
		   )
		 ),
		"Hands are weird");



class FilterEventsWithBigNegDeviation : public G3Module { 
 public:
 FilterEventsWithBigNegDeviation(std::string frb_event_key, 
				 std::string out_frb_event_key  ) :
  frb_event_key_(frb_event_key), out_frb_event_key_(out_frb_event_key)
  {}

  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string out_frb_event_key_;
};

EXPORT_G3MODULE("frbutils", 
		FilterEventsWithBigNegDeviation, 
		(init<std::string, std::string>(
		(bp::arg("frb_event_key"),
		 bp::arg("out_frb_event_key")
		 ))),
		"Hands are weird");






class FilterForDetectorPair : public G3Module {
 public:
  FilterForDetectorPair(std::string frb_event_key,
			std::string frb_event_out_key);
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string frb_event_out_key_;  
};
EXPORT_G3MODULE("frbutils", 
		FilterForDetectorPair, 
		(
		 bp::init<std::string, std::string>
		 ( (
		    bp::arg("frb_event_key"),
		    bp::arg("frb_event_out_key")
		    )
		   )
		 ),
		"Hands are weird");



class FilterNumDets : public G3Module {
 public:
  FilterNumDets(std::string frb_event_key,
		std::string frb_event_out_key,
		size_t n_dets);
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string frb_event_key_;
  std::string frb_event_out_key_;  
  size_t n_dets_;
};
EXPORT_G3MODULE("frbutils", 
		FilterNumDets, 
		(
		 bp::init<std::string, std::string, size_t>
		 ( ( bp::arg("frb_event_key"),
		     bp::arg("frb_event_out_key"),
		     bp::arg("n_dets")
		     )
		   )
		 ),
		"Hands are weird");



class SmartAmpFilter : public G3Module {
 public:
  SmartAmpFilter(std::string frb_event_key,
		std::string frb_event_out_key,
		std::string variance_path, 
		 double n_sigma_away,
		 double max_amp_err,
		 bool invert_result);
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
	std::string frb_event_key_;
	std::string frb_event_out_key_;  
	std::string variance_path_;
	double n_sigma_away_;
	double max_amp_err_;
	bool invert_result_;
};

EXPORT_G3MODULE("frbutils", 
		SmartAmpFilter, 
		(bp::init<std::string, std::string, std::string, double, double, bool>
		 ( ( bp::arg("frb_event_key"),
		     bp::arg("frb_event_out_key"),
		     bp::arg("variance_path"),
		     bp::arg("n_sigma_away"),
		     bp::arg("max_amp_err"),
		     bp::arg("invert_result")
			 )
			 )
			),
		"Hands are weird");


class MakeOtherThresholds : public G3Module {
public:
	MakeOtherThresholds(std::string frb_event_key,
			    std::string frb_event_out_key,
			    double event_threshold,
			    double other_threshold
		) : frb_event_key_(frb_event_key),
		    frb_event_out_key_(frb_event_out_key),
		    event_threshold_(event_threshold),
		    other_threshold_(other_threshold) {}
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	std::string frb_event_key_;
	std::string frb_event_out_key_;  
	double event_threshold_;
	double other_threshold_;
};

EXPORT_G3MODULE("frbutils", 
		MakeOtherThresholds, 
		(
		 bp::init<std::string, std::string, double, double>
		 ( ( bp::arg("frb_event_key"),
		     bp::arg("frb_event_out_key"),
		     bp::arg("event_threshold"),
		     bp::arg("other_threshold")
		     )
		   )
		 ),
		"Hands are weird");
