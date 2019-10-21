#pragma once

#include <string>
#include <vector>
#include <G3Timestream.h>
#include <coordinateutils/G3SkyMap.h>

void bin_tod_pol(const G3Timestream & ts, double pol_angle, double pol_eff,
			 const std::vector<int> & map_indices, double det_wgt,
			 G3SkyMap & T, G3SkyMap & Q, G3SkyMap & U);

void bin_w(const G3Timestream & ts, const std::vector<int> & map_indices, double det_wgt,
	   G3SkyMapWeights & W);

void bin_w_pol(const G3Timestream & ts, double pol_angle, double pol_eff,
	       const std::vector<int> & map_indices, double det_wgt,
	       G3SkyMapWeights & W);

void fill_tod_pol(const std::vector<int> & map_indices, 
		  const G3SkyMap & T, const G3SkyMap & Q, const G3SkyMap & U,
		  double pol_angle, double pol_eff,
		  G3Timestream & ts);

void bin_tod_pol_rot(const G3Timestream & ts, double pol_angle, double pol_eff,
		     const std::vector<int> & map_indices, double det_wgt,
		     const std::vector<double> & pol_rot,
		     G3SkyMap & T, G3SkyMap & Q, G3SkyMap & U);

void bin_w_pol_rot(const G3Timestream & ts, double pol_angle, double pol_eff,
		   const std::vector<int> & map_indices, double det_wgt,
		   const std::vector<double> & pol_rot,
		   G3SkyMapWeights & W);

void fill_tod_pol_rot(const std::vector<int> & map_indices, 
		      const G3SkyMap & T, const G3SkyMap & Q, const G3SkyMap & U,
		      const std::vector<double> & pol_rot,
		      double pol_angle, double pol_eff,
		      G3Timestream & ts);

void fill_tod_pol_interp_2d(const std::vector<double> & alphas,
			    const std::vector<double> & deltas,
			    const G3SkyMap & T, const G3SkyMap & Q,
			    const G3SkyMap & U,
			    double pol_angle, double pol_eff,
			    G3Timestream & ts);

void fill_tod_pol_rot_interp_2d(const std::vector<double> & alphas,
				const std::vector<double> & deltas,
				const G3SkyMap & T, const G3SkyMap & Q,
				const G3SkyMap & U,
				const std::vector<double> & pol_rot,
				double pol_angle, double pol_eff,
				G3Timestream & ts);
