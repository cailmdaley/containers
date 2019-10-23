#include <pybindings.h>
#include <G3Logging.h>

#include <string>
#include <vector>

#include <mapmaker/polarizationutils.h>

void bin_tod_pol(const G3Timestream & ts, double pol_angle, double pol_eff,
		 const std::vector<int> & map_indices, double det_wgt,
		 G3SkyMap & T, G3SkyMap & Q, G3SkyMap & U) {
	if (det_wgt == 0) {
		return;
	}

	StokesVector sc;
	MuellerMatrix pcm;
	set_stokes_coupling(pol_angle, pol_eff, 1.0, sc);
	for (size_t i=0; i<ts.size(); i++){
		if (ts[i] == 0)
			continue;
		T[map_indices.at(i)] += det_wgt * sc.t * ts[i];
		Q[map_indices.at(i)] += det_wgt * sc.q * ts[i];
		U[map_indices.at(i)] += det_wgt * sc.u * ts[i];
	}
}

void bin_w(const G3Timestream & ts, const std::vector<int> & map_indices, double det_wgt,
	       G3SkyMapWeights & W){
	if (det_wgt == 0) {
		return;
	}
	for (size_t i=0; i<ts.size(); i++){
		(*W.TT)[map_indices.at(i)] += det_wgt;
	}
}

void bin_w_pol(const G3Timestream & ts, double pol_angle, double pol_eff,
	       const std::vector<int> & map_indices, double det_wgt,
	       G3SkyMapWeights & W){
	if (det_wgt == 0) {
		return;
	}
	StokesVector sc;
	MuellerMatrix pcm;
	set_stokes_coupling(pol_angle, pol_eff, 1.0, sc);
	fill_mueller_matrix_from_stokes_coupling(sc, pcm);
	pcm *= det_wgt;
	for (size_t i=0; i<ts.size(); i++){
		W[map_indices.at(i)] += pcm;
	}
}

void fill_tod_pol(const std::vector<int> & map_indices, 
		  const G3SkyMap & T, const G3SkyMap & Q, const G3SkyMap & U,
		  double pol_angle, double pol_eff,
		  G3Timestream & ts) {
	StokesVector sc;
	set_stokes_coupling(pol_angle, pol_eff, 1.0, sc);
	for (size_t i=0; i<ts.size(); i++){
		ts[i] = T.at(map_indices.at(i)) * sc.t +
			Q.at(map_indices.at(i)) * sc.q +
			U.at(map_indices.at(i)) * sc.u;
	}
	ts.units = T.units;
}

void bin_tod_pol_rot(const G3Timestream & ts, double pol_angle, double pol_eff,
		     const std::vector<int> & map_indices, double det_wgt,
		     const std::vector<double> & pol_rot,
		     G3SkyMap & T, G3SkyMap & Q, G3SkyMap & U){
	StokesVector sc;
	for (size_t i=0; i<ts.size(); i++){
		if (ts[i] == 0)
			continue;
		set_stokes_coupling(pol_angle+pol_rot[i], pol_eff, 1.0, sc);
		T[map_indices.at(i)] += det_wgt * sc.t * ts[i];
		Q[map_indices.at(i)] += det_wgt * sc.q * ts[i];
		U[map_indices.at(i)] += det_wgt * sc.u * ts[i];
	}
}

void bin_w_pol_rot(const G3Timestream & ts, double pol_angle, double pol_eff,
		   const std::vector<int> & map_indices, double det_wgt,
		   const std::vector<double> & pol_rot,
		   G3SkyMapWeights & W){
	StokesVector sc;
	MuellerMatrix pcm;

	for (size_t i=0; i<ts.size(); i++){
		set_stokes_coupling(pol_angle+pol_rot[i], pol_eff, 1.0, sc);
		fill_mueller_matrix_from_stokes_coupling(sc, pcm);
		pcm *= det_wgt;
		W[map_indices.at(i)] += pcm;
	}
}

void fill_tod_pol_rot(const std::vector<int> & map_indices, 
		      const G3SkyMap & T, const G3SkyMap & Q, const G3SkyMap & U,
		      const std::vector<double> & pol_rot,
		      double pol_angle, double pol_eff,
		      G3Timestream & ts) {
	StokesVector sc;
	for (size_t i=0; i<ts.size(); i++){
		set_stokes_coupling(pol_angle+pol_rot[i], pol_eff, 1.0, sc);
		ts[i] = T.at(map_indices.at(i)) * sc.t +
			Q.at(map_indices.at(i)) * sc.q +
			U.at(map_indices.at(i)) * sc.u;
	}
	ts.units = T.units;
}

void fill_tod_pol_interp_2d(const std::vector<double> & alphas,
			    const std::vector<double> & deltas,
			    const G3SkyMap & T, const G3SkyMap & Q,
			    const G3SkyMap & U,
			    double pol_angle, double pol_eff,
			    G3Timestream & ts) {
	StokesVector sc;
	set_stokes_coupling(pol_angle, pol_eff, 1.0, sc);
	for (size_t i=0; i<ts.size(); i++){
		std::vector<long> pixels;
		std::vector<double> weights;
		T.get_interp_pixels_weights(alphas[i], deltas[i], pixels, weights);
		ts[i] = T.get_interp_precalc(pixels, weights) * sc.t +
			Q.get_interp_precalc(pixels, weights) * sc.q +
			U.get_interp_precalc(pixels, weights) * sc.u;
	}
	ts.units = T.units;
}

void fill_tod_pol_rot_interp_2d(const std::vector<double> & alphas,
				const std::vector<double> & deltas,
				const G3SkyMap & T, const G3SkyMap & Q,
				const G3SkyMap & U,
				const std::vector<double> & pol_rot,
				double pol_angle, double pol_eff,
				G3Timestream & ts) {
	StokesVector sc;
	for (size_t i=0; i<ts.size(); i++){
		set_stokes_coupling(pol_angle+pol_rot[i], pol_eff, 1.0, sc);
		std::vector<long> pixels;
		std::vector<double> weights;
		T.get_interp_pixels_weights(alphas[i], deltas[i], pixels, weights);
		ts[i] = T.get_interp_precalc(pixels, weights) * sc.t +
			Q.get_interp_precalc(pixels, weights) * sc.q +
			U.get_interp_precalc(pixels, weights) * sc.u;
	}
	ts.units = T.units;
}
