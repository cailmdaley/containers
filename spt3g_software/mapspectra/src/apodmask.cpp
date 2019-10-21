#include <pybindings.h>
#include <coordinateutils/FlatSkyMap.h>

/*
 * Zero all pixels with other zero pixels within n_clean pixels.
 */

static FlatSkyMapPtr
erode(FlatSkyMapPtr in_m, int n_clean, int bthresh)
{
	FlatSkyMapPtr m(new FlatSkyMap(*in_m));
	for (size_t d = 0; d < n_clean; d++){
		// Iteratively zero pixels with too few non-zero pixels
		// bordering it. This zeroes out to a radius equal to the
		// number of iterations.
		
		FlatSkyMap n(*m); // Grab current state
		for (size_t y = 1; y < m->shape()[1] - 1; y++) {
			for (size_t x = 1; x < m->shape()[0] - 1; x++) {
				// Count non-zero pixels neighboring this one
				int border = 0;
				for (int x_s = -1; x_s <= 1; x_s++) {
					for (int y_s = -1; y_s <= 1; y_s++) {
						if (y_s == 0 && x_s == 0)
							continue;
						if (n(x+x_s, y+y_s) != 0)
							border++;
					}
				}
				if (border < bthresh)
					(*m)(x,y) = 0;
			}
		}
	}
	return m;
}

static FlatSkyMapPtr
erode_cardinal(FlatSkyMapPtr in_m, int n_clean, int bthresh)
{
	FlatSkyMapPtr m(new FlatSkyMap(*in_m));
	for (size_t d = 0; d < n_clean; d++){
		// Iteratively zero pixels with too few non-zero pixels
		// bordering it. This zeroes out to a radius equal to the
		// number of iterations.
		
		const FlatSkyMap n(*m); // Grab current state
		for (size_t y = 1; y < m->shape()[1] - 1; y++) {
			for (size_t x = 1; x < m->shape()[0] - 1; x++) {
				// Skip pixels already zero
				if (((const FlatSkyMap &)(*m))(x,y) == 0)
					continue;

				// Count non-zero pixels neighboring this one
				int border = 0;
				if (n(x+1, y) != 0)
					border++;
				if (n(x-1, y) != 0)
					border++;
				if (n(x, y+1) != 0)
					border++;
				if (n(x, y-1) != 0)
					border++;
				if (border < bthresh)
					(*m)(x,y) = 0;
			}
		}
	}
	return m;
}

PYBINDINGS("mapspectra") {
	namespace bp = boost::python;
	bp::def("erode", erode, (bp::arg("r_clean"), bp::arg("bthresh")=8),
	    "Iteratively zero all pixels within r_clean pixels of a zero "
	    "pixel. The optional parameter bthresh sets the number of allowed "
	    "non-zero neighboring pixels at each step. The default (8) "
	    "requires that all neighboring pixels be non-zero.");
	bp::def("erode_cardinal", erode_cardinal, (bp::arg("r_clean"), bp::arg("bthresh")=4),
	    "Iteratively zero all pixels within r_clean pixels of a zero "
	    "pixel. The optional parameter bthresh sets the number of allowed "
	    "non-zero neighboring pixels at each step. The default (4) "
	    "requires that all neighboring pixels be non-zero.  This differs from erode"
	    " in that it does not include the diagonal pixels.");
}

