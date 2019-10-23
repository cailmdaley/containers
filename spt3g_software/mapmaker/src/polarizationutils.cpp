#include <G3Units.h>

#include <mapmaker/polarizationutils.h>


void fill_mueller_matrix_from_stokes_coupling( const StokesVector & stokes_coupling, 
					       MuellerMatrix & pcm){
	pcm.tt = stokes_coupling.t * stokes_coupling.t;
	pcm.tq = stokes_coupling.t * stokes_coupling.q;
	pcm.tu = stokes_coupling.t * stokes_coupling.u;
	pcm.qq = stokes_coupling.q * stokes_coupling.q;
	pcm.qu = stokes_coupling.q * stokes_coupling.u;
	pcm.uu = stokes_coupling.u * stokes_coupling.u;
}

void
set_stokes_coupling(double pol_ang, double pol_eff, double cal_constant,
                    StokesVector & stokes_coupling)
{
	stokes_coupling.t =  cal_constant;
	stokes_coupling.q =  cos( pol_ang/G3Units::rad * 2. ) * pol_eff * cal_constant/(2.0-pol_eff);
	stokes_coupling.u =  sin( pol_ang/G3Units::rad * 2. ) * pol_eff * cal_constant/(2.0-pol_eff);

	stokes_coupling.q = fabs(stokes_coupling.q) > 1e-12 ? stokes_coupling.q : 0.0;
	stokes_coupling.u = fabs(stokes_coupling.u) > 1e-12 ? stokes_coupling.u : 0.0;
}
