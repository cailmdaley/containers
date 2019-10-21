#pragma once

#include <vector>
#include <G3Frame.h>
#include <G3Module.h>
#include <G3Timestream.h>
#include <G3Map.h>


void make_low_pass_filter(double sample_rate,  double low_pass_cutoff, size_t len,
			  std::vector<std::complex< double > > & filter);

void fft_filter_mem_friendly(const G3TimestreamMap & time_data, 
			     const G3VectorComplexDouble & filter_data,
			     G3TimestreamMap & out_data,
			     
			     bool decon_timestreams = false,
			     G3MapDoubleConstPtr time_constants = G3MapDoubleConstPtr() 
	);

