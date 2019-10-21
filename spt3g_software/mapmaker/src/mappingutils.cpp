#include <pybindings.h>
#include <serialization.h>

#include <sstream>
#include <G3Logging.h>

#include <string>
#include <vector>

#include <mapmaker/mappingutils.h>
#include <mapmaker/polarizationutils.h>


#define COS cos
#define SIN sin

std::vector<double> get_unit_xyz(double ra, double dec) {
	double x = COS(ra / G3Units::rad) *  COS(dec/G3Units::rad);
	double y = SIN(ra / G3Units::rad) *  COS(dec/G3Units::rad);
	double z = SIN(dec / G3Units::rad);
	return {x, y, z};
}

void get_unit_xyz_vec(const std::vector<double> & ra_lst,
		      const std::vector<double> & dec_lst,
		      std::vector<double> & x,
		      std::vector<double> & y,
		      std::vector<double> & z)
{
	x = std::vector<double>(ra_lst.size());
	y = std::vector<double>(ra_lst.size());
	z = std::vector<double>(ra_lst.size());
	for (unsigned int i=0; i < ra_lst.size(); i++){
		std::vector<double> xyz = get_unit_xyz(ra_lst[i], dec_lst[i]);
		x[i] = xyz[0];
		y[i] = xyz[1];
		z[i] = xyz[2];
	}
}

//assumes out_mask has been allocated
void make_point_source_mask(const std::vector<double> & ra_lst,
			    const std::vector<double> & dec_lst,
			    const std::vector<double> & radius_lst,
			    bool mask_out_of_bounds_pixels,
			    G3SkyMap & out_mask)
{
	std::vector<double> ptsrc_x;
	std::vector<double> ptsrc_y;
	std::vector<double> ptsrc_z;

	g3_assert(ra_lst.size() == dec_lst.size());
	g3_assert(ra_lst.size() == radius_lst.size());

	get_unit_xyz_vec(ra_lst, dec_lst, ptsrc_x, ptsrc_y, ptsrc_z);

	std::vector<double> cos_rad_dist = std::vector<double>(radius_lst.size());
	for (size_t i=0; i < radius_lst.size(); i++){
		cos_rad_dist[i] = COS(radius_lst[i]/G3Units::rad);
	}

	for (size_t j = 0; j < out_mask.size(); j++) {
		std::vector<double> radec = out_mask.pixel_to_angle(j);
		std::vector<double> vec = get_unit_xyz(radec[0], radec[1]);
		double v = 0;
		for (unsigned int i=0; i < ptsrc_x.size(); i++)
			v +=  ! ! ((ptsrc_x[i] * vec[0] + ptsrc_y[i] * vec[1] + ptsrc_z[i] * vec[2]) > cos_rad_dist[i]);
		if (!!v)
			out_mask[j] = 1.0;
	}

	out_mask.overflow = mask_out_of_bounds_pixels;
}

void make_determinant_map(G3SkyMapWeightsConstPtr W, G3SkyMapPtr & Dout) {
	Dout = (*W).TT->Clone(true);
	for (size_t i = 0; i < W->TT->size(); i++) {
		double det = W->at(i).det();
		if (det != 0)
			(*Dout)[i] = det;
	}
}

void remove_weight(G3SkyMapConstPtr T, G3SkyMapConstPtr Q,
		   G3SkyMapConstPtr U, G3SkyMapWeightsConstPtr W,
		   G3SkyMapPtr & Tout, G3SkyMapPtr & Qout, G3SkyMapPtr & Uout) {
	g3_assert(T->size() == Q->size());
	g3_assert(T->size() == U->size());
	g3_assert(W->TT->size() == T->size());
	if (! T->is_weighted )log_warn("T map may have already had the weight removed");
	if (! Q->is_weighted )log_warn("Q map may have already had the weight removed");
	if (! U->is_weighted )log_warn("U map may have already had the weight removed");
	Tout = T->Clone(true);
	Qout = Q->Clone(true);
	Uout = U->Clone(true);
	for (size_t i = 0; i < T->size(); i++) {
		StokesVector v((*Tout)[i], (*Qout)[i], (*Uout)[i]);
		v /= W->at(i);
	}
	Tout->is_weighted = false;
	Qout->is_weighted = false;
	Uout->is_weighted = false;
}

static boost::python::object
remove_weight_py(G3SkyMapConstPtr T, G3SkyMapConstPtr Q, G3SkyMapConstPtr U,
    G3SkyMapWeightsConstPtr W)
{
	G3SkyMapPtr outT, outQ, outU;
	remove_weight(T, Q, U, W, outT, outQ, outU);
	
	return boost::python::make_tuple(outT, outQ, outU);
}


void remove_weight_t(G3SkyMapConstPtr T, 
		     G3SkyMapWeightsConstPtr W,
		     G3SkyMapPtr & Tout) {
	g3_assert(W->TT->size() == T->size());
	if (! T->is_weighted )log_warn("T map may have already had the weight removed");
	Tout = T->Clone(true);
	*Tout /= *W->TT;
	Tout->is_weighted = false;
}

namespace bp = boost::python;

void mappingutils_pybindings(void){
	bp::def("make_point_source_mask_cpp", make_point_source_mask,
		(bp::arg("ra_lst"),
		 bp::arg("dec_lst"),
		 bp::arg("radius_lst"),
		 bp::arg("mask_out_of_bounds_pixels"),
		 bp::arg("map")
			)
		);

	bp::def("remove_weight", remove_weight_py,
		(bp::arg("T"),
		 bp::arg("Q"),
		 bp::arg("U"),
		 bp::arg("W")
		),
                "Remove weight from a set of maps, returning unweighted T,Q,U");
	
	bp::def("remove_weight_t_cpp", remove_weight_t, "Removes the weight for an unpolarized map while properly accounting for division by 0.");
	bp::def("make_determinant_map_cpp", make_determinant_map);
}

