#pragma once

#include <string>
#include <deque>
#include <G3Frame.h>
#include <G3Module.h>
#include <G3Map.h>
#include <G3Vector.h>
#include <G3Timestream.h>
#include <G3Data.h>
#include <G3Logging.h>

#include <todfilter/polyutils.h>
#include <todfilter/fftutils.h>

FilterMaskConstPtr get_scan_filter_mask(G3FramePtr scan_frame,
					std::string pointing_key, 
					boost::shared_ptr<const G3SkyMap > point_source_mask );

#define FILTER_MASK_INJECTOR_DOCSTR \
"Converts a point source mask in map space into a mask in timestream space for filtering.\n\n" \
"\n" \
"pntsrc_mask [G3SkyMapConstPtr]: The point source mask.  1 means masked, 0 means not masked.\n\n" \
"filter_mask_key: The path where the filter mask should be stored.\n\n"

class FilterMaskInjector : public G3Module{
 public:
	FilterMaskInjector(std::string pntsrc_mask_id, 
			   std::string filter_mask_key, 
			   std::string pointing_key);
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
	
private:
	std::string pointing_key_;
	std::string mask_key_;
	std::string pntsrc_mask_id_; 
	G3SkyMapConstPtr pntsrc_mask_; 
	SET_LOGGER("FilterMaskInjector");
};



#define FFT_FILTER_DOCSTR \
"Runs the FFT based filter with the specified filter.\n\n" \
"\n\n" \
"ts_in_key [->G3TimestreamMap]  input timestreams\n\n" \
"ts_out_key [->G3TimestreamMap] output timestreams\n\n" \
"fft_filter_key [->G3MapVectorComplexDouble]  path to the specified filter\n\n"

class FftFilter : public G3Module{
public:
	FftFilter(std::string ts_in_key, std::string ts_out_key, 
		  std::string fft_filter_key, std::string padding_key) : 
		ts_in_key_(ts_in_key), ts_out_key_(ts_out_key), fft_filter_key_(fft_filter_key),
		padding_key_(padding_key){}
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
	//lowpass filter
	//detector deconvolution
private:
	std::string ts_in_key_;
	std::string ts_out_key_;
	std::string fft_filter_key_;
	std::string padding_key_;
	SET_LOGGER("FftFilter");
};


#define ROLLING_MEAN_FILTER_DOCSTR \
	"ts_in_key [->G3TimestreamMap]  input timestreams\n\n" \
	"ts_out_key [->G3TimestreamMap]  output timestreams\n\n" \
	"filter_width:  the width of the filter to use with mean filtering"


class RollingMeanFilter : public G3Module{
 public:
  RollingMeanFilter(std::string ts_in_key, std::string ts_out_key, 
	    int filter_width) : 
  ts_in_key_(ts_in_key), ts_out_key_(ts_out_key), filter_width_(filter_width){}
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
  //lowpass filter
  //detector deconvolution
 private:
  std::string ts_in_key_;
  std::string ts_out_key_;
  int filter_width_;
  SET_LOGGER("RollingMeanFilter");
};




#define MASKED_POLY_HPF_DOCSTR \
"Does a masked poly and high pass filter with a linear least squares fitter.\n\n" \
"\n\n" \
"in_ts_map_key: path to input timestream map\n\n" \
"out_ts_map_key: path to output timestream map\n\n" \
"poly_order [int]: order of polynomial to filter.  poly_order < 0 means no poly filter.\n\n"\
"high_pass_freq_cutoff [double]: frequency to cutoff the HPF. If negative, no mHPF is performed.\n\n"\
"is_masked [bool] whether or not it is masked\n\n" \
"mask_key : path to a filter mask, not used if is_masked is false\n\n"

class MaskedPolyHpf : public G3Module{
public:
	MaskedPolyHpf(std::string in_ts_map_key, std::string out_ts_map_key,
		      int poly_order, double high_pass_freq_cutoff, 
		      bool is_masked, std::string mask_key, 
		      std::string sample_rate_override_key);
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	bool is_masked_;
	int poly_order_;
	double high_pass_freq_cutoff_;
	std::string in_ts_map_key_;
	std::string out_ts_map_key_;
	std::string mask_key_;
	
	std::string sample_rate_override_key_;

	SET_LOGGER("MaskedPolyHpf");
};
G3_POINTERS(MaskedPolyHpf);
G3_SERIALIZABLE(MaskedPolyHpf, 1);


#define MASKED_NOTCH_FILTER_DOCSTR \
"Does a masked notch filter by taking out    \n" \
"sines and cosines at specified frequencies. \n" \
"This is not FFT-based notching but          \n" \
"masked-high-pass-based notching.            \n" \
"\n\n\n" \
"input_tsm_key  (str): key name of the input timestream map.  \n\n" \
"output_tsm_key (str): key name of the output timestream map. \n\n" \
"frequencies [G3VectorDouble]: list of frequencies to use for notching. \n\n" \
"is_masked [bool]: whether or not timestreams are masked. \n\n" \
"mask_key (str)  : key name of the filter mask.           \n\n"

class MaskedNotchFilter: public G3Module{
public:
    MaskedNotchFilter(
        std::string input_tsm_key,
        std::string output_tsm_key,
        const std::vector<double> & frequencies,
        bool is_masked,
        std::string mask_key);
    void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
    std::string input_tsm_key_;
    std::string output_tsm_key_;
    const std::vector<double> frequencies_;
    bool is_masked_;
    std::string mask_key_;
    SET_LOGGER("MaskedNotchFilter");
};
G3_POINTERS(MaskedNotchFilter);
G3_SERIALIZABLE(MaskedNotchFilter, 1);
