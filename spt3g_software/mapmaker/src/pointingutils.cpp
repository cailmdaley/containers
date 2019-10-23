#include <pybindings.h>

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>

#include <mapmaker/pointingutils.h>

#include <iostream>

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

#define COS cos
#define SIN sin
#define ASIN asin
#define ATAN2 atan2

using namespace G3Units;


// az_out and el out are a 1d represenatation of 2d information where bolo index is the slow changing, time index is fast

void bs_pointing_to_bolo_delta_alpha(const std::vector<angle_t> & alpha_bs,
				     const std::vector<angle_t> & delta_bs,
				     const BolometerPropertiesMap & bolo_props,
				     MapCoordReference coord_sys,
				     G3MapVectorDouble & alpha_out,
				     G3MapVectorDouble & delta_out,
				     const G3VectorString & dets
	){
	
	//log_warn("bs_pointing_to_bolo_delta_alpha is deprecated");
	size_t n_times = alpha_bs.size();
	size_t n_bolos = bolo_props.size();
	size_t i;

	//boresight parameters
	std::vector<double> sin_alpha_bs(n_times);
	std::vector<double> cos_alpha_bs(n_times);
	std::vector<double> sin_delta_bs(n_times);
	std::vector<double> cos_delta_bs(n_times);

	//detector parameters
	std::vector<double> sin_alpha_in(n_bolos);
	std::vector<double> cos_alpha_in(n_bolos);
	std::vector<double> sin_delta_in(n_bolos);
	std::vector<double> cos_delta_in(n_bolos);
	std::vector<double> tan_delta_in(n_bolos);
	//right now the map maker is not set up to handle rotatations of the focal plane relative
	//to the coordinate system used.  it could be set up to do this, but I don't think we need it
	//g3_assert(coord_sys == Equatorial || coord_sys == Local);
	double delta_mul = coord_sys == Local ? 1.0f : -1.0f;

	g3_assert(alpha_bs.size() == delta_bs.size());

	//ra dec adjusted pointing
	#ifdef OPENMP_FOUND
	#pragma omp parallel for
	#endif
	for (i=0; i < n_times; i++){
		sin_alpha_bs[i] = SIN(alpha_bs[i]/rad);
		cos_alpha_bs[i] = COS(alpha_bs[i]/rad);
		sin_delta_bs[i] = SIN(delta_bs[i]/rad);
		cos_delta_bs[i] = COS(delta_bs[i]/rad);
	}

	// OpenMP can't see through std::map iterators, so build a list first
	/**
	std::vector<std::string> dets;
	for (auto it=bolo_props.begin(); it != bolo_props.end(); it++)
	  dets.push_back(it->first);
	**/
	#ifdef OPENMP_FOUND
	#pragma omp parallel for
	#endif
	for (i=0; i < dets.size(); i++){
		sin_alpha_in[i] = SIN(bolo_props.at(dets[i]).x_offset/rad);
		cos_alpha_in[i] = COS(bolo_props.at(dets[i]).x_offset/rad);
		sin_delta_in[i] = delta_mul * SIN(bolo_props.at(dets[i]).y_offset/rad);
		cos_delta_in[i] = COS(bolo_props.at(dets[i]).y_offset/rad);
		tan_delta_in[i] = sin_delta_in[i]/cos_delta_in[i];
	}

	for (size_t i=0; i < dets.size(); i++){
		alpha_out[dets[i]] = G3VectorDouble(n_times, 0);
		delta_out[dets[i]] = G3VectorDouble(n_times, 0);
	}
	// loop over the bolometers

	#ifdef OPENMP_FOUND
	#pragma omp parallel for
	#endif
	for(i=0; i < dets.size(); i++){
		std::vector<double> & alpha_ts = alpha_out[dets[i]];
		std::vector<double> & delta_ts = delta_out[dets[i]];
		for (size_t j=0; j< n_times; j++){
			delta_ts[j] = ASIN(sin_delta_bs[j]*cos_delta_in[i]*cos_alpha_in[i]
			                   + cos_delta_bs[j]*sin_delta_in[i]  ) * rad;
			alpha_ts[j] = -ATAN2( (- sin_alpha_bs[j]*cos_delta_bs[j]*cos_alpha_in[i] 
					       - cos_alpha_bs[j]*sin_alpha_in[i] 
					       + sin_alpha_bs[j]*sin_delta_bs[j]*tan_delta_in[i] ),
			                      ( cos_alpha_bs[j]*cos_delta_bs[j]*cos_alpha_in[i] 
						- sin_alpha_bs[j]*sin_alpha_in[i]  
						- cos_alpha_bs[j]*sin_delta_bs[j]*tan_delta_in[i] )
			                    ) * rad;
		}
	}
}

void bs_pointing_to_bolo_delta_alpha_one_samp(double alpha_bs, double delta_bs,
					      const BolometerPropertiesMap & bolo_props,
					      const G3VectorString & dets,
					      G3MapDouble & alpha_out,
					      G3MapDouble & delta_out) {
	std::vector<angle_t> alpha_bs_vec(1, alpha_bs);
	std::vector<angle_t> delta_bs_vec(1, delta_bs);


	G3MapVectorDouble alpha_out_vec;
	G3MapVectorDouble delta_out_vec;
	bs_pointing_to_bolo_delta_alpha(alpha_bs_vec, delta_bs_vec, bolo_props, 
					MapCoordReference::Equatorial,
					alpha_out_vec, delta_out_vec, dets);

	for (auto it = alpha_out_vec.begin(); it != alpha_out_vec.end(); it++) {
		alpha_out[it->first] = alpha_out_vec[it->first][0];
		delta_out[it->first] = delta_out_vec[it->first][0];
	}
}


namespace bp = boost::python;
//PYBINDINGS("mapmaker"){
void pointingutils_pybindings(void){
	bp::def("bs_pointing_to_bolo_delta_alpha", bs_pointing_to_bolo_delta_alpha,
		(bp::arg("alpha_bs"), bp::arg("delta_bs"),
		 bp::arg("bolo_props"), bp::arg("coord_sys"),
		 bp::arg("alpha_out"), bp::arg("delta_out"),
		 bp::arg("dets")
			),
		"Converts the boresight pointing to alpha/delta coordinates. Works only for local and equatorial coordinates and relies on the fact that we have a telescope very close to the south pole when doing the az/el offsets.  If dets is non zero length, only calculates for the dets in dets."
		);
	bp::def("bs_pointing_to_bolo_delta_alpha_one_samp", 
		bs_pointing_to_bolo_delta_alpha_one_samp);
}
