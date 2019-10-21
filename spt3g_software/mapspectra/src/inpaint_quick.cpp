#include <pybindings.h>

#include <coordinateutils/FlatSkyMap.h>

#include <mapspectra/inpaint_quick.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>


int inpaint_map_laplace(FlatSkyMapConstPtr mask, int n_iters,
			FlatSkyMapPtr map){
	if (!mask->IsCompatible(*map)) {
		log_fatal("Map and mask not same dimensions");
	}

	int i,j;
	size_t nx = mask->shape()[0];
	size_t ny = mask->shape()[1];
	std::vector<int> pixel_x_inds;
	std::vector<int> pixel_y_inds;
	
	for (i=0; i<ny; i++){
		for (j=0; j<nx; j++){
			if (mask->at(j, i)){
				pixel_x_inds.push_back(j);
				pixel_y_inds.push_back(i);
				(*map)(j, i) = 0;
			}
		}
	}
        size_t n_pixels = pixel_x_inds.size();
        std::vector<double> tmp_pixel_vals(n_pixels);
	
	//actually do the inpainting
	for (i=0; i<n_iters; i++){
		//if (i %100 == 0) printf("On loop %d\n", i);
		for (j=0; j<n_pixels; j++){
			int i_x = pixel_x_inds[j];
			int i_y = pixel_y_inds[j];
			
			double tmp_sum = 0;
			double nsum = 0;
			
			//edge adjacent
			if (i_x - 1 >= 0){
				tmp_sum += map->at(i_x - 1, i_y);
				nsum += 1.0;
			}
			if (i_x + 1 < nx){
				tmp_sum += map->at(i_x + 1, i_y);
				nsum += 1.0;
			}
			if (i_y - 1 >= 0){
				tmp_sum += map->at(i_x, i_y - 1);
				nsum += 1.0;
			}
			if (i_y + 1 < ny){
				tmp_sum += map->at(i_x, i_y + 1);
				nsum += 1.0;
			}
			tmp_sum /= nsum;
			tmp_pixel_vals[j] = tmp_sum;
		}
		for (j=0; j<n_pixels; j++){
			(*map)(pixel_x_inds[j], pixel_y_inds[j]) = tmp_pixel_vals[j];
		}
	}

	return 0;
}

namespace bp = boost::python;
PYBINDINGS("mapspectra"){
	bp::def("inpaint_map_laplace", inpaint_map_laplace, 
		(bp::arg("mask"), bp::arg("n_iters"), bp::arg("map")),
		"Inpaints over the masked region.  For the mask 1 means that it is inpainted over.  0 means nothing happens to that pixel."
		);
}
