#include <todfilter/polyutils.h>
#include <G3Timestream.h>
#include <G3Map.h>
#include <iostream>
#include <string>
//poly_and_mhpf_filter_ts_data

#define TS_LEN 2000
int main() {

	G3TimestreamMapPtr tsm(new G3TimestreamMap);
	for (size_t i=0; i < 1000; i++){
		std::string s = std::to_string(i);
		(*tsm)[s]= G3TimestreamPtr(new G3Timestream(TS_LEN, 1));
	}
	FilterMaskPtr fm = make_empty_filter_mask(tsm);


	/**
	for (auto it = fm->has_masked_pixels.begin(); it != fm->has_masked_pixels.end(); it++){
		it->second = 1;
		fm->pixel_mask[it->first] = G3VectorInt(TS_LEN,0);
	}
	**/
	for (size_t i=0; i < 30; i++){
		G3TimestreamMapPtr tso = poly_and_mhpf_filter_ts_data(tsm, fm, 1, 0.001, 10);
	}
	return 0;
}
