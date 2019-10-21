#include <pybindings.h>
#include <serialization.h>
#include <todfilter/numericutils.h>
#include <G3Map.h>

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fftw3.h>

#include <iostream>
#include <vector>
#include <complex>

#ifdef OPENMP_FOUND
#include <omp.h>
#define FFTW_OMP
#endif




float hann_window(int index, int n_indices){
	g3_assert(index >= 0 && index < n_indices);
	float s =  sinf(M_PI * index / (n_indices - 1.0f));
	return s * s;
}

float hamming_window(int index, int n_indices){
	g3_assert(index >= 0 && index < n_indices);
	return 0.54 - .46 * cosf( 2.0 * M_PI * index / (n_indices - 1.0f));
}


std::vector<double> fft_freqs(double sample_frequency, int fftlen){
	std::vector<double> freqs(fftlen);
	if (fftlen % 2 == 0) {
		for (int i=0; i < fftlen/2; i++) freqs[i] = i * sample_frequency / fftlen;
		for (int i = fftlen/2; i < fftlen; i++) freqs[i] = ( i - fftlen ) * sample_frequency / fftlen;
	} else {
		for (int i=0; i < fftlen/2+1; i++) freqs[i] = i * sample_frequency / fftlen;
		for (int i = fftlen/2+1; i < fftlen; i++) freqs[i] = ( i - 2.0 * (fftlen/2) -1 ) * sample_frequency / fftlen;
	}
	return freqs;
}

//sptsz uses       lpfilt[i]=exp(-1*pow(freqs[i]/lowpass,6));


inline double filt_func(double f, double low_pass_cutoff){
	return exp( -1.0 * pow(f / low_pass_cutoff, 12) );
}


void make_low_pass_filter(double sample_rate,  double low_pass_cutoff, size_t len,
			  std::vector<std::complex< double > > & filter){
	//
	//auto filt_func = [&] (double f) -> std::complex<double> { return exp( -1.0 * pow(f / low_pass_cutoff, 6) )  ;};
	filter= std::vector<std::complex< double > >(len);
	std::vector<double> freqs = fft_freqs(sample_rate, len);
	for (size_t i=0; i < len; i++){
		filter[i] = filt_func(freqs[i], low_pass_cutoff);
	}
}




void bulk_fft_forward(const G3TimestreamMap & time_data,
		      size_t filter_len,  G3MapVectorComplexDouble & freq_data){
	size_t ts_len = time_data.begin()->second->size();
	g3_assert(ts_len <= filter_len);
	size_t filter_delta = filter_len - ts_len;
	std::vector<std::string> id_lst(time_data.size());
	size_t i=0;

	size_t complex_filter_len = filter_len/2 + 1;

	for (auto it=time_data.begin(); it != time_data.end(); it++, i++){
		id_lst[i] = it->first;
	}
	freq_data = G3MapVectorComplexDouble();
	for (auto iter = time_data.begin(); iter != time_data.end(); iter++)
		freq_data[iter->first] = std::vector< std::complex< double >  >(complex_filter_len, 0);
#ifdef FFTW_OMP
#pragma omp parallel
#endif
	{
		double * timed;
		fftw_complex * freqd;
		fftw_plan p;
#ifdef FFTW_OMP
#pragma omp critical(fftw)
#endif
		{
			timed = (double*) fftw_malloc(sizeof(double) * filter_len);
			freqd = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * complex_filter_len);
			p = fftw_plan_dft_r2c_1d(filter_len, timed, freqd, 
						 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
		}
#ifdef FFTW_OMP
#pragma omp for 
#endif
		for (i=0; i < id_lst.size(); i++){
			memcpy(timed, &( (*time_data.at(id_lst[i]))[0]), ts_len * sizeof(double) );
			memset(timed+ts_len, 0, sizeof(double) * filter_delta);
			fftw_execute(p);
			memcpy((void *)&(freq_data[id_lst[i]][0]), freqd, complex_filter_len * sizeof(fftw_complex ) );
		}
#ifdef FFTW_OMP
#pragma omp critical(fftw)
#endif
		{ 
			fftw_destroy_plan(p);
			fftw_free(timed); 
			fftw_free(freqd);
		}
	}
}



G3MapVectorDoublePtr get_psd(G3TimestreamMapPtr time_data_in, size_t padding_amount,
			     G3VectorDoublePtr window_func, double norm_fac) {

	G3TimestreamMapPtr time_data(new G3TimestreamMap);
	for (auto iter = time_data_in->begin(); iter != time_data_in->end(); iter++){
		(*time_data)[iter->first] = G3TimestreamPtr(new G3Timestream(*(iter->second)));
	}

	size_t ts_len = time_data->begin()->second->size();
	
	for (auto it = time_data->begin(); it != time_data->end(); it++){
		G3Timestream &ts = *(it->second);

		double ts_mean=0;
		for (size_t i=0; i < ts_len; i++) {
			ts_mean += ts[i];
		}
		ts_mean= ts_mean/ts_len;

		for (size_t i=0; i < ts_len; i++) {
			ts[i] =(ts[i]-ts_mean) *(*window_func)[i];

		}
	}

	double mean_window_func_val = 0;	
	for (size_t i=0; i < ts_len; i++) {
		mean_window_func_val += (*window_func)[i] * (*window_func)[i];
	}
	mean_window_func_val /= ts_len;

	G3Map< std::string, std::vector< std::complex< double > > > freq_data;
	bulk_fft_forward(*time_data, padding_amount, freq_data);
	
	size_t freq_len = freq_data.begin()->second.size();
	
	G3MapVectorDoublePtr out_data(new G3MapVectorDouble());
	for (auto it = freq_data.begin(); it != freq_data.end(); it++){
		std::vector< std::complex< double > > & ft = it->second;
		G3VectorDouble psd(freq_len);
		for (size_t i=0; i < freq_len; i++) {
			double a = std::abs(ft[i]);
			psd[i] = 2.0 * a * a * norm_fac / ( mean_window_func_val * ts_len);
		}
		//corrects for the sections that aren't repeated
		psd[0] /= 2.0;

		if (padding_amount % 2 == 0 ) 
			psd[ freq_len - 1 ] /= 2.0;

		(*out_data)[it->first] = psd;
	}
	return out_data;
}


/**
   f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
   f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
 **/


/**  
     If time_data and out_data point to the same object memory is not allocated.

**/

//uses too much memory for the complex filter
//also weirdly slow.  Fuck
void fft_filter_mem_friendly(const G3TimestreamMap & time_data, 
			     const G3VectorComplexDouble & filter_data,
			     G3TimestreamMap & out_data,
			     bool decon_timestreams,
			     G3MapDoubleConstPtr time_constants 
	){
	g3_assert(time_data.size() > 0);
	double sample_rate = time_data.begin()->second->GetSampleRate();
	
	int filter_len = filter_data.size();
	int ts_len = time_data.begin()->second->size();
	G3Map< std::string, std::vector< std::complex< double > > > freq_data;
	std::vector<double> hamming_window_vec(ts_len);  
	const double i_filt_len_d = 1.0/((double)filter_len);
	int complex_filter_len = filter_data.size() / 2 + 1;	


	if (&out_data != &time_data) {
		out_data = G3TimestreamMap();
		for (auto i = time_data.begin(); i != time_data.end(); i++)
			out_data[i->first] = G3TimestreamPtr(new G3Timestream(*i->second));
	}
	else {
		log_fatal("Attempting to do fft filter in place.  That's not how we do things here");
	}
	
	//creates the hamming window
	for (int i=0; i < ts_len; i++) hamming_window_vec[i] = hamming_window(i, ts_len);  
	
	//creates the list of ids because iterators don't play nice with openmp
	std::vector<std::string> id_lst(out_data.size());
	size_t i=0;
	for (auto it=out_data.begin(); it != out_data.end(); it++, i++){
		id_lst[i] = it->first;
	}
	
	size_t filter_delta = filter_len - ts_len;
	
	// if we are deconvolving precalculate the 
	std::vector<double> freqs;
	if (decon_timestreams){
		g3_assert(time_constants);
		freqs = fft_freqs(sample_rate, filter_len);
	}
	
	
#ifdef FFTW_OMP
#pragma omp parallel
#endif
	{
		
		//allocates the fft plans
		double * timed;
		fftw_complex * freqd;
		fftw_plan forward_p, backward_p;
		
#ifdef FFTW_OMP
#pragma omp critical(fftw)
#endif
		{    
			timed = (double*) fftw_malloc(sizeof(double) * filter_len);
			freqd = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 
							    complex_filter_len);
			forward_p  = fftw_plan_dft_r2c_1d(filter_len, timed, freqd,
							  FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
							  //FFTW_MEASURE | FFTW_DESTROY_INPUT);
			backward_p = fftw_plan_dft_c2r_1d(filter_len, freqd, timed, 
							  FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
							  //FFTW_MEASURE | FFTW_DESTROY_INPUT);
		}
		
		
		//and now the fun
#ifdef FFTW_OMP
#pragma omp for 
#endif
		for (i=0; i < id_lst.size(); i++) {
			G3TimestreamPtr ts = out_data[id_lst[i]];
			//scale by hamming window
			for (int j=0; j < ts_len; j++) (*ts)[j] *= hamming_window_vec[j];
			// fft forward
			memcpy(timed, &((*ts)[0]), ts_len * sizeof(double) );
			memset(timed+ts_len, 0, sizeof(double) * filter_delta);
			fftw_execute(forward_p);
			
	  //apply filter
			std::complex< double > * ffted_ts = (std::complex< double > *) freqd;
			for (int j=0; j < complex_filter_len; j++) 
				ffted_ts[j] *= filter_data[j];
			
			//probably apply timestream decon
			if (decon_timestreams) {
				std::complex<double> jmaginary_number(0.0,1.0);
				double tc = time_constants->at(id_lst[i]);
				for (int j=0; j < complex_filter_len; j++) {
					ffted_ts[j] *= ( 1.0 +  2.0 * M_PI * jmaginary_number * freqs[j] * tc );
				}
			}
			//invert
			fftw_execute(backward_p);
			//invert window
			for (int j=0; j < ts_len; j++) 
				timed[j] *= i_filt_len_d/hamming_window_vec[j];
			memcpy(&((*ts)[0]), timed, ts_len * sizeof(double));
		}
		
#ifdef FFTW_OMP
#pragma omp critical(fftw)
#endif
	{
		fftw_destroy_plan(forward_p);
		fftw_destroy_plan(backward_p);
		fftw_free(timed); 
		fftw_free(freqd);
	}
	
	}
}



void fft_filter_mem_friendly_w_multi_filters(
    G3TimestreamMap           &  timestream_map,
    G3VectorDouble            &  window_function,
    G3MapVectorComplexDouble  &  filter_map,
    G3TimestreamMap           &  filtered_ts_map)
    
    // This function is a modified version of
    // the function fft_filter_mem_friendly above.
    // The main difference is that this function
    // can apply different Fourier space filters to different timestreams,
    // as opposed to applying the same filter to all timestreams.

{
    double sample_rate = timestream_map.begin()->second->GetSampleRate();
    int    ts_len      = timestream_map.begin()->second->size();
    int    filter_len  = filter_map.begin()->second.size();
    
    int          complex_filter_len = filter_len / 2 + 1;
    const double i_filt_len_d       = 1.0 / ((double)filter_len);
    
    
    if (&filtered_ts_map != &timestream_map)
    {
        filtered_ts_map = G3TimestreamMap();
        for (auto iter=timestream_map.begin();
                  iter!=timestream_map.end(); iter++)
        {
            filtered_ts_map[iter->first] = 
                G3TimestreamPtr(new G3Timestream(*iter->second));
        }
    }
    else
    {
        log_fatal("Attempting to do fft filter in place. Not allowed!");
    }
    
    
    std::vector<std::string> bolonames(filtered_ts_map.size());
    size_t i = 0;
    for (auto iter=filtered_ts_map.begin();
              iter!=filtered_ts_map.end(); iter++, i++)
    {
        bolonames[i] = iter->first;
    }
    
    
    size_t filter_delta = filter_len - ts_len;
    
    
#ifdef FFTW_OMP
#pragma omp parallel
#endif
    
    {
        double      * timed;
        fftw_complex* freqd;
        fftw_plan     forward_p, backward_p;
        
#ifdef FFTW_OMP
#pragma omp critical(fftw)
#endif
    
        {
            timed = (double*)       fftw_malloc(sizeof(double) *
                                                filter_len);
            freqd = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *
                                                complex_filter_len);
            
            forward_p  = fftw_plan_dft_r2c_1d(
                             filter_len, timed, freqd,
                             FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            backward_p = fftw_plan_dft_c2r_1d(
                             filter_len, freqd, timed,
                             FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
        }
        
        
#ifdef FFTW_OMP
#pragma omp for
#endif
        
        for (i=0; i<bolonames.size(); i++)
        {
            G3TimestreamPtr ts = filtered_ts_map[bolonames[i]];
            
            for (int j=0; j<ts_len; j++)
            {
                (*ts)[j] *= window_function[j];
            }
            
            memcpy(timed, &((*ts)[0]), ts_len * sizeof(double));
            memset(timed+ts_len, 0, sizeof(double) * filter_delta);
            
            fftw_execute(forward_p);
            std::complex< double > * ffted_ts =
                (std::complex< double > *) freqd;
            
            for (int j=0; j<complex_filter_len; j++)
            {
                ffted_ts[j] *= filter_map[bolonames[i]][j];
            }
            
            fftw_execute(backward_p);
            
            for (int j=0; j<ts_len; j++)
            {
                timed[j] *= i_filt_len_d / window_function[j];
            }
            
            memcpy(&((*ts)[0]), timed, ts_len * sizeof(double));
        }
        
        
#ifdef FFTW_OMP
#pragma omp critical(fftw)
#endif
        
        fftw_destroy_plan(forward_p);
        fftw_destroy_plan(backward_p);
        fftw_free(timed); 
        fftw_free(freqd);
    
    }

}




namespace bp = boost::python;
PYBINDINGS("todfilter"){
	bp::def("fft_filter_mem_friendly", 
		fft_filter_mem_friendly,
		( bp::arg("time_data"),
		  bp::arg("filter_data"),
		  bp::arg("out_data"),
		  bp::arg("decon_timestreams"),
		  bp::arg("time_constants")
			),
		"This is meant to be a memory-friendly version of the FFT filter.  "
		"Applies a complex frequency space filter to the time data.  "
		"Right now there is stub code for deconvolving time constant effects with the system.  "
		"The code has several inefficiencies that could be improved.  "
		"The timestream deconvolution also hasn't been tested.  "
		"This will probably need to be updated for a DAN circuit model.  "
		"By default it applies a hamming window to the data since there are no points where "
		"that window function is zero.  It removes the window after filtering.  "
		"\n\nThe filter is defined by the complex values of filter_data.  Where each index "
		"in filter data corresponds to the frequency of the DFT'd data. numpy's fftfreqs should "
		"give you that.  If the filter is longer than the input time data, it zero pads the "
		"time data by the appropriate amount."
		);

	bp::def("fft_filter_mem_friendly_w_multi_filters", 
		fft_filter_mem_friendly_w_multi_filters,
		( bp::arg("timestream_map"),
		  bp::arg("window_function"),
		  bp::arg("filter_map"),
		  bp::arg("filtered_ts_map")
			),
                "This funciton is a modified version of "
                "fft_filter_mem_friendly. "
                "Instead of applying one filter to "
                "all the timestreams in a timestream map, "
                "this function allows one to apply different filters "
                "to different timestreams. "
		);

	bp::def("get_psd_pybinding_grr", get_psd, 
		(bp::arg("time_data"),
		 bp::arg("padding_amount"),
		 bp::arg("window_func"),
		 bp::arg("norm_fac")
			), 
		"This function is parseval's theorem obeying.\n"
		);

	bp::def("fft_freqs", fft_freqs);

}

