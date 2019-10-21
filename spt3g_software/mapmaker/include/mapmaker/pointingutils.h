#pragma once
#include <vector>

#include <G3Frame.h>
#include <G3Module.h>
#include <G3Timestream.h>
#include <G3Map.h>
#include <G3Vector.h>

#include <calibration/BoloProperties.h>
#include <coordinateutils/G3SkyMap.h>

typedef double angle_t;

void bs_pointing_to_bolo_delta_alpha(const std::vector<double> & alpha_in,
				     const std::vector<double> & delta_in,
				     const BolometerPropertiesMap & bolo_props,
				     MapCoordReference coord_sys,
				     G3MapVectorDouble & alpha_out, 
				     G3MapVectorDouble & delta_out,
				     const G3VectorString & dets);

void pointingutils_pybindings(void);
