#pragma once
#include <pybindings.h>
#include <G3Frame.h>
#include <G3Map.h>
#include <G3Timestream.h>
#include <vector>
#include <string>
#include <stdlib.h>
#include <G3Logging.h>

#include <coordinateutils/G3SkyMap.h>

#include <boost/python.hpp>

void make_point_source_mask(const std::vector<double> & ra_lst,
                            const std::vector<double> & dec_lst,
                            const std::vector<double> & radius_lst,
                            bool mask_out_of_bounds_pixels,
                            G3SkyMap & map);

void remove_weight(G3SkyMapConstPtr T, G3SkyMapConstPtr Q,
		   G3SkyMapConstPtr U, G3SkyMapWeightsConstPtr W,
		   G3SkyMapPtr & Tout, G3SkyMapPtr & Qout, G3SkyMapPtr & Uout);

void mappingutils_pybindings();
