#pragma once
#include <coordinateutils/G3SkyMap.h>

void set_stokes_coupling(double pol_ang, double pol_eff, double cal_constant,
			 StokesVector & stokes_coupling);
void fill_mueller_matrix_from_stokes_coupling( const StokesVector & stokes_coupling,
					       MuellerMatrix & pcm);
