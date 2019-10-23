import numpy as np
import camb
from scipy import special

pars = camb.read_ini('/home/tcrawfor/.local/lib/python3.6/site-packages/camb/planck_2018.ini')
tfs = camb.get_transfer_functions(pars)
temp = tfs.get_cmb_transfer_data()
tlq = temp.get_transfer()
ell = tlq[0]
k = tlq[1]
transfer_t = tlq[2]
pk = 2.100549e-9*(k/0.05)**(0.9660499-1.)
clcmb = np.zeros(len(ell))
for i in np.arange(len(ell)):
    clcmb[i] = 4.*np.pi*np.trapz(1./k*transfer_t[i,:]**2*pk,k)
clmut = np.zeros(len(ell))
#mu_ave = 2.3e-8
mu_ave = 2e-8
qd = 0.084
for i in np.arange(len(ell)):
    clmut[i] = 24./5./np.pi*mu_ave*np.trapz(1./k*2.*np.pi**2*pk*transfer_t[i,:]*np.exp(-(k/qd)**2)*special.spherical_jn(ell[i],k*14000.*0.994),k)

ellcmb = np.arange(np.max(ell)+1)
clcmb_int = np.interp(ellcmb,ell,clcmb)
clmut_int = np.interp(ellcmb,ell,clmut)

clmm_n_1ukarcmin = np.zeros(len(ellcmb))+(1e-6/2.73*np.pi/180./60.)**2

f00_v_lmax_1uk = np.zeros(len(ellcmb))
f00_v_lmax_1uk[0:3] = 1e-12
for i in np.arange(len(ellcmb)):
    if i > 2:
        f00_v_lmax_1uk[i] = np.sum((2.*ellcmb[2:i+1]+1.)*clmut_int[2:i+1]**2/(clmm_n_1ukarcmin[2:i+1]*clcmb_int[2:i+1]))

# try to reproduce PIXIE curve from Cabass et al. as cross-check
clmm_n_pixie = np.zeros(len(ellcmb)) + 1e12
clmm_n_pixie[0:500] = 4.*np.pi*1e-16*np.exp((ellcmb[0:500]/84)**2)
f00_v_lmax_pixie = np.zeros(len(ellcmb))
f00_v_lmax_pixie[0:3] = 1e-12
for i in np.arange(len(ellcmb)):
    if i > 2:
        f00_v_lmax_pixie[i] = np.sum((2.*ellcmb[2:i+1]+1.)*clmut_int[2:i+1]**2/(clmm_n_pixie[2:i+1]*clcmb_int[2:i+1]))

# CMB-S4 baseline (from eyeballing ILC)
clmm_n_cmbs4 = np.zeros(len(ellcmb))+(8e-6/2.73*np.pi/180./60.)**2
f00_v_lmax_cmbs4 = np.zeros(len(ellcmb))
f00_v_lmax_cmbs4[0:3] = 1e-12
for i in np.arange(len(ellcmb)):
    if i > 2:
        f00_v_lmax_cmbs4[i] = np.sum((2.*ellcmb[2:i+1]+1.)*clmut_int[2:i+1]**2/(clmm_n_cmbs4[2:i+1]*clcmb_int[2:i+1]))





