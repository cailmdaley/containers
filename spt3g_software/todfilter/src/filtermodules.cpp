#include <pybindings.h>
#include <todfilter/filtermodules.h>



FilterMaskConstPtr get_scan_filter_mask(G3FramePtr scan_frame,
					std::string pointing_key, 
					G3SkyMapConstPtr point_source_mask ){
	G3MapVectorIntConstPtr pointing_indices = scan_frame->Get<G3MapVectorInt>(pointing_key);
	FilterMaskConstPtr filter_mask = make_filter_mask(point_source_mask, pointing_indices);
	return filter_mask;
}



void FftFilter::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	log_debug("Entering Process");
	if (frame->type == G3Frame::Scan){
		g3_assert(frame->Has(ts_in_key_));
		g3_assert(frame->Has(fft_filter_key_));
		G3TimestreamMapPtr out_data = boost::make_shared<G3TimestreamMap>();
		G3TimestreamMapConstPtr time_data = frame->Get<G3TimestreamMap>(ts_in_key_);
		G3VectorComplexDoubleConstPtr filter_data = frame->Get<G3VectorComplexDouble>(fft_filter_key_);
		
		int ts_len = time_data->begin()->second->size();
		int filter_len = ts_len;
		if (!padding_key_.empty()){
			filter_len = frame->Get<G3Int>(padding_key_)->value;
		}
		g3_assert(time_data);
		g3_assert(out_data);
		g3_assert(filter_data);
		fft_filter_mem_friendly(*time_data, *filter_data, *out_data, false, G3MapDoubleConstPtr());
		frame->Put(ts_out_key_, out_data);
	}
	out.push_back(frame);
}



//want to change it to accept a string that's the map id
//


FilterMaskInjector::FilterMaskInjector(std::string pntsrc_mask_id,
				       std::string filter_mask_key, 
				       std::string pointing_key):
	pointing_key_(pointing_key), mask_key_(filter_mask_key), 
	pntsrc_mask_id_(pntsrc_mask_id){}


void FilterMaskInjector::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	log_debug("Entering Process");
	if (frame->type == G3Frame::Scan){
		if (! pntsrc_mask_) {
			log_fatal("Point source mask has not been put through yet.");
		}
		if (! frame->Has(mask_key_)){
			FilterMaskConstPtr fmaskptr = get_scan_filter_mask(frame, pointing_key_, 
									   pntsrc_mask_);
			frame->Put(mask_key_, fmaskptr);
		}
	} else if (frame->type == G3Frame::Map) {
		if (frame->Get<G3String>("Id")->value == pntsrc_mask_id_)
			pntsrc_mask_ = frame->Get<G3SkyMap>("T");
	}
	out.push_back(frame);
}




MaskedPolyHpf::MaskedPolyHpf(std::string in_ts_map_key,
    std::string out_ts_map_key, int poly_order, double high_pass_freq_cutoff, 
			     bool is_masked, std::string mask_key, 
			     std::string sample_rate_override_key) :
	is_masked_(is_masked), poly_order_(poly_order), 
	high_pass_freq_cutoff_(high_pass_freq_cutoff),
	in_ts_map_key_(in_ts_map_key), out_ts_map_key_(out_ts_map_key), mask_key_(mask_key),
	sample_rate_override_key_(sample_rate_override_key){}


void MaskedPolyHpf::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	log_debug("Entering Process");
	if ( frame->type == G3Frame::Scan){
		//do scan things
		g3_assert(frame->Has(in_ts_map_key_));
		G3TimestreamMapConstPtr in_ts = frame->Get<G3TimestreamMap>(in_ts_map_key_);
		G3TimestreamMapPtr out_ts;// = boost::make_shared<G3TimestreamMap>();
		double sample_frequency = in_ts->GetSampleRate();

		if (sample_rate_override_key_.size() != 0) {
			sample_frequency = 
				(frame->Get<G3Double>(sample_rate_override_key_))->value;
		}

		FilterMaskConstPtr filter_mask;
		if (is_masked_){
			log_debug("Grabing mask\n");
			filter_mask = frame->Get<FilterMask>(mask_key_);
		}  else {
			filter_mask = make_empty_filter_mask(in_ts);
		}
		out_ts = poly_and_mhpf_filter_ts_data(in_ts, filter_mask, 
						      sample_frequency, high_pass_freq_cutoff_, 
						      poly_order_);
		frame->Put(out_ts_map_key_, out_ts);
	}
	out.push_back(frame);
}



MaskedNotchFilter::MaskedNotchFilter(
    std::string input_tsm_key,
    std::string output_tsm_key,
    const std::vector<double> & frequencies,
    bool is_masked,
    std::string mask_key):

    input_tsm_key_(input_tsm_key),
    output_tsm_key_(output_tsm_key),
    frequencies_(frequencies),
    is_masked_(is_masked),
    mask_key_(mask_key){}

void MaskedNotchFilter::Process(
    G3FramePtr frame,
    std::deque<G3FramePtr> &out){
    
    if (frame->type == G3Frame::Scan){
        g3_assert(frame->Has(input_tsm_key_));
        G3TimestreamMapConstPtr input_tsm = 
            frame->Get<G3TimestreamMap>(input_tsm_key_);
        G3TimestreamMapPtr output_tsm;
        
        FilterMaskConstPtr filter_mask;
        if (is_masked_){
            filter_mask = frame->Get<FilterMask>(mask_key_);
        } else {
            filter_mask = make_empty_filter_mask(input_tsm);
        }
        output_tsm = notch_filter_lls(input_tsm, filter_mask, frequencies_);
        frame->Put(output_tsm_key_, output_tsm);
    }
    
    out.push_back(frame);
}
    



void RollingMeanFilter::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	if ( frame->type == G3Frame::Scan){
		G3TimestreamMapConstPtr in_ts = frame->Get<G3TimestreamMap>(ts_in_key_);
		G3TimestreamMapPtr out_ts = boost::make_shared<G3TimestreamMap>();
		rolling_mean_filter_ts_data(in_ts, filter_width_, out_ts);
		frame->Put(ts_out_key_, out_ts);
	}
	out.push_back(frame);
}




namespace bp = boost::python;

PYBINDINGS("todfilter"){
	/**
  bp::class_<FftFilter, bp::bases<G3Module>,
	     boost::shared_ptr<FftFilter>,
	     boost::noncopyable >("FftFilter",
				  bp::init<std::string, std::string, std::string, std::string> 
				  ( FFT_FILTER_DOCSTR,
				    (bp::arg("ts_in_key"), 
				     bp::arg("ts_out_key"), 
				     bp::arg("fft_filter_key"),
				     bp::arg("padding_key")
					    )) 
		     )
    .def("Process", &FftFilter::Process)
    ;
  bp::implicitly_convertible< boost::shared_ptr<FftFilter>, G3ModulePtr>();	
	**/


  bp::class_<RollingMeanFilter, bp::bases<G3Module>,
	     boost::shared_ptr<RollingMeanFilter>,
	     boost::noncopyable >("RollingMeanFilter",
				  bp::init<std::string, std::string, int> 
				  ( ROLLING_MEAN_FILTER_DOCSTR, (bp::arg("ts_in_key"), 
				     bp::arg("ts_out_key"), 
				     bp::arg("filter_width"))) 
		     )
    .def("Process", &RollingMeanFilter::Process)
    ;
  bp::implicitly_convertible< boost::shared_ptr<RollingMeanFilter>, G3ModulePtr>();	


  bp::class_<MaskedPolyHpf, bp::bases<G3Module>,
	     boost::shared_ptr<MaskedPolyHpf>,
	     boost::noncopyable >("MaskedPolyHpf",
				  bp::init<std::string, std::string, int, double, bool, 
				  std::string, std::string>
				  ( MASKED_POLY_HPF_DOCSTR,
				    (bp::arg("in_ts_map_key"), bp::arg("out_ts_map_key"), 
				    bp::arg("poly_order"), bp::arg("high_pass_freq_cutoff") = 0,
				     bp::arg("is_masked") = false, bp::arg("mask_key") = "",
				     bp::arg("sample_rate_override_key") = ""
				    ))
				  )
    .def("Process", &MaskedPolyHpf::Process)
    ;
  bp::implicitly_convertible< boost::shared_ptr<MaskedPolyHpf>, G3ModulePtr>();	


    bp::class_<MaskedNotchFilter, bp::bases<G3Module>,
        boost::shared_ptr<MaskedNotchFilter>,
        boost::noncopyable>("MaskedNotchFilter",
                            bp::init<std::string,
                                     std::string,
                                     const std::vector<double>,
                                     bool,
                                     std::string>
                            (MASKED_NOTCH_FILTER_DOCSTR,
                             (bp::arg("input_tsm_key"),
                              bp::arg("output_tsm_key"), 
                              bp::arg("frequencies"),
                              bp::arg("is_masked")=false,
                              bp::arg("mask_key")="")))
    .def("Process", &MaskedNotchFilter::Process)
    ;
    bp::implicitly_convertible<boost::shared_ptr<MaskedNotchFilter>, G3ModulePtr>();

                   
  bp::class_<FilterMaskInjector, bp::bases<G3Module>,
	     boost::shared_ptr<FilterMaskInjector>,
	     boost::noncopyable >("FilterMaskInjector", 
				  bp::init<std::string, std::string, std::string>
				  (FILTER_MASK_INJECTOR_DOCSTR,
				   ( bp::arg("point_src_mask_id"), 
				     bp::arg("filter_mask_key"), 
				     bp::arg("pointing_key")
				      )) 
				  )
    .def("Process", &FilterMaskInjector::Process)
    ;
  bp::implicitly_convertible< boost::shared_ptr<FilterMaskInjector>, G3ModulePtr>();	

}
