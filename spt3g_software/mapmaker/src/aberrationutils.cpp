#include <pybindings.h>
#include <mapmaker/aberrationutils.h>

#include <G3Units.h>
#include <coordinateutils/G3SkyMap.h>

#include <boost/geometry.hpp>

using namespace G3Units;
namespace bg = boost::geometry;


void aberrate_delta_alpha(G3MapVectorDouble alpha_in,
                          G3MapVectorDouble delta_in,
                          MapCoordReference coord_sys,
                          G3MapVectorDouble & alpha_out,
                          G3MapVectorDouble & delta_out,
                          double vel, double alpha, double delta
                          ){
	g3_assert(coord_sys == Equatorial || coord_sys == Local);
	g3_assert(alpha_in.size() == delta_in.size());

	double beta = vel/(300000.0*km/s);
	double gamma = 1.0/sqrt(1.0+pow(beta,2));
	double sin_delta = sin(delta/rad);
	double cos_delta = cos(delta/rad);

	int n_times = alpha_in.begin()->second.size();
	std::vector<double> sin_delta_in(n_times);
	std::vector<double> cos_delta_in(n_times);
	std::vector<double> cos_alpha_diff(n_times);

	for (auto it=alpha_in.begin(); it != alpha_in.end(); it++){
		alpha_out[it->first] = G3VectorDouble(n_times, 0);
		delta_out[it->first] = G3VectorDouble(n_times, 0);
	}

	// loop over the bolometers
	for (auto it=alpha_in.begin(); it != alpha_in.end(); it++){
		std::vector<double> & alpha_ts = alpha_out[it->first];
		std::vector<double> & delta_ts = delta_out[it->first];
		const std::vector<double> & alpha_old_ts = alpha_in.at(it->first);
		const std::vector<double> & delta_old_ts = delta_in.at(it->first);
		for (size_t j=0; j< n_times; j++){
			cos_alpha_diff[j] = cos((alpha-alpha_old_ts[j])/rad);
			sin_delta_in[j] = sin(delta_old_ts[j]/rad);
			cos_delta_in[j] = cos(delta_old_ts[j]/rad);
			// Planck 2013 XXVII Eq. 2
			double correction_mag = (gamma-1.0)*(sin_delta*sin_delta_in[j]*cos_alpha_diff[j]
			                         +cos_delta*cos_delta_in[j])+gamma*beta;
			bg::model::point<double, 3, bg::cs::spherical<bg::radian> >
				dir(alpha_old_ts[j]/rad, delta_old_ts[j]/rad+M_PI/2.0, 1.0);
			bg::model::point<double, 3, bg::cs::spherical<bg::radian> >
				velocity(alpha, delta+M_PI/2.0, correction_mag);
			bg::model::point<double, 3, bg::cs::cartesian > dir_cart;
			bg::model::point<double, 3, bg::cs::cartesian > vel_cart;
			bg::transform(dir, dir_cart);
			bg::transform(velocity, vel_cart);
			bg::add_point(dir_cart, vel_cart);
			bg::transform(dir_cart, dir);

			alpha_ts[j] = bg::get<0>(dir)*rad;
			delta_ts[j] = bg::get<1>(dir)*rad-M_PI/2.0;
		}
	}
}

void aberrationutils_pybindings(void){
  namespace bp = boost::python;
	bp::def("aberrate_delta_alpha", aberrate_delta_alpha,
		(bp::arg("alpha_in"),
		 bp::arg("delta_in"),
		 bp::arg("coord_sys"),
		 bp::arg("alpha_out"),
		 bp::arg("delta_out"),
		 bp::arg("vel"),
		 bp::arg("alpha"),
		 bp::arg("delta")),
		"Applies relativistic aberration to detector pointing"
		);
}
