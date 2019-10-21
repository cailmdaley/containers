import os, sys, scipy, hashlib, glob, subprocess, imp, pdb
import pickle as pk
import healpy as hp
import numpy as np
import pylab as pl
import datetime
import glob
from spt3g import core
from spt3g.lensing.map_spec_utils import MapSpectraTEB
from spt3g import lensing as sl

# ----------------------
# Cuts
# ----------------------
lmax_sinv = 4000  # cutoff the filtered signal at this value

# When evalutating phi, cutoff the integral between 'lmax_cinv' and 'lx_cut'
lmax_cinv = 3000  # cutoff the cinv-filtered fields at this value
lx_cut = 450  # mask modes below this value in lx
nsim = 500

# T->P leakage
tq_leak = 0.0050
tu_leak = -0.0083
# tq_leak= -0.0050
# tu_leak= 0.0083
# polarization calibration.  Note, this is applied in addition to Tcal.
# tcal = 1.0 ## use default value, since we take care of Tcal in the beam
pcal = 1.048

# ----------------------
# Names
# ----------------------
ivfs_prefix = "run08"
qest_prefix = "mf100_lx450_lmax3000"

# Specify which sims are used for qecl's
mc_sims = np.arange(0, nsim)
mc_sims_mf = mc_sims[0:100]
mc_sims_var = mc_sims[100:500]
mc_sims_unl = mc_sims_var
mc_sims_qcr_mc = mc_sims_var
mc_sims_n0 = mc_sims
mc_sims_n1 = mc_sims[0:50]


# cmbs
core.log_info("---- get cmbs", unit="ParameterFile")
bdir = "/spt/user/panz/lens100d_new_pipeline/"
if not os.path.exists(bdir):
    core.log_notice(
        "bdir does not exist. The code assumes you are running it on the grid."
    , unit="ParameterFile")
    bdir = (
        glob.glob(
            os.path.abspath(os.path.dirname(__file__)).split("code")[0]
            + "lensing_*"
        )[0].split(".tar.gz")[0]
        + "/"
    )

scan = imp.load_source("scan", bdir + "scripts/make_cmbs_scan.py")
cmbs = scan.cmbs

# maps (13 deg, 2 arcmin, 390 pixels)
reso_rad = np.double(scan.reso / 60 * np.pi / 180.0)  # 2 arcmin
npix = int(scan.nx)  # 390
lmax_theoryspec = scan.lmax
# 6000.  Cut off the theory spectrum at this value

cl_unl = sl.map_spec_utils.get_camb_scalcl(
    prefix="planck_wp_highL", lmax=lmax_theoryspec
)
cl_unl["BB"] = np.zeros(lmax_theoryspec + 1)
cl_len = sl.map_spec_utils.get_camb_lensedcl(
    prefix="planck_wp_highL", lmax=lmax_theoryspec
)

core.log_info("---- get beams", unit="ParameterFile")
blpk150 = pk.load(
    open(
        bdir + "inputs/beam/composite_beam_2012_13_bb_runlist_v3_150.pkl", "rb"
    ),
    encoding="latin1",
)
map_cal_factor = blpk150["cal"]  # composite calibration number, = 0.8926

bl150 = blpk150["B_ell_hybrid"] / map_cal_factor
tf_beam_cal = (
    MapSpectraTEB(
        ffts=3 * [np.ones([npix, npix // 2 + 1])],
        map_nx=npix,
        map_ny=npix,
        dx=reso_rad,
        dy=reso_rad,
    )
    * bl150
).get_l_masked(lmax=lmax_theoryspec)


core.log_info("---- get TF", unit="ParameterFile")
tf_file = (
    bdir
    + "inputs/tf2d/tf2d_sptpol_poly4_lowpass4000_nobeam_ngrid390_reso2am_proj5_1000iter.npy"
)
tf150 = (
    MapSpectraTEB(
        ffts=3 * [np.load(tf_file) * tf_beam_cal.get_pixel_window()],
        map_nx=npix,
        map_ny=npix,
        dx=reso_rad,
        dy=reso_rad,
    )
    * bl150
).get_l_masked(lmax=lmax_theoryspec)

core.log_info("---- get fgndlib", unit="ParameterFile")
fgndlib150 = sl.lens100d_interface.fgnd.ForegroundT(
    npix,
    npix,
    reso_rad,
    reso_rad,
    lmax_theoryspec,
    bdir + "qest/fgnds",
    asrc=10.0,
    acib=5.0,
    atsz=5.0,
    pcib=0.8,
    seed=44,
)  ##
halfs_dir = bdir + "data/map/20140418_ra23h30dec55_lens100d/halfs/150/"

core.log_info("---- get half_noise", unit="ParameterFile")
obs150_len = sl.obs.ObsHalfNoise(
    tf_beam_cal,
    sl.obs.SumTQU([scan.cmblib_len_t1p1_scan, fgndlib150]),
    halfs_dir,
    tq=tq_leak,
    tu=tu_leak,
    thresh=1000.0 * core.G3Units.uK,
    pcal=pcal,
)
obs150_unl = sl.obs.ObsHalfNoise(
    tf_beam_cal,
    sl.obs.SumTQU([scan.cmblib_unl_t1p1_scan, fgndlib150]),
    halfs_dir,
    tq=tq_leak,
    tu=tu_leak,
    thresh=1000.0 * core.G3Units.uK,
    pcal=pcal,
)
# Note, no foregrounds needed in obs150_len_t2,
# since it is only used for estimating the N1 bias
obs150_len_t2_nofg = sl.obs.ObsHomogeneousNoise(
    tf_beam_cal,
    sl.obs.SumTQU([scan.cmblib_len_t2p1_scan]),
    thresh=1000.0 * core.G3Units.uK,
)
# for calculating the N1 bias as: len_t2 - len-nofg
obs150_len_nofg = sl.obs.ObsHomogeneousNoise(
    tf_beam_cal,
    sl.obs.SumTQU([scan.cmblib_len_t1p1_scan]),
    thresh=1000.0 * core.G3Units.uK,
)
obslibs = [obs150_len, obs150_len_nofg, obs150_len_t2_nofg]

# ivfs
core.log_info("---- ivfs", unit="ParameterFile")
# used for plotting purposes only
apod150 = np.load(bdir + "inputs/mask/apod_surv_five_150.npy")
mask150 = (
    np.load(bdir + "inputs/mask/mask_surv_five_150.npy")
    * np.load(bdir + "inputs/mask/mask_clust10.npy")
    * np.load(bdir + "inputs/mask/mask_srccr10.npy")
)  ##

core.log_info("---- ivfs:ninv150", unit="ParameterFile")
ninv150 = sl.map_spec_utils.make_tqumap_wt(
    map_nx=npix,
    map_ny=npix,
    dx=reso_rad,
    dy=reso_rad,
    ninv=obs150_len.get_ninv(),
    ninv_dcut=1.0e-5,
    nlev_tp=(7.0, 7.0),
    mask=mask150,
)
# map-space noise
ninvfilt150 = sl.cinv_utils.opfilt_teb.NoiseInverseFilter(tf150, ninv150)

core.log_info("---- ivfs:sinvfilt150", unit="ParameterFile")
cl_len_filt = {}
for k in cl_len.keys():
    cl_len_filt[k] = cl_len[k][0 : lmax_sinv + 1]
del cl_len_filt["TE"]  # theoretical spectrum
# includes signal and fourier-space noise
fgnd_cl = fgndlib150.get_clfg(lmax=lmax_sinv)
total_cl = {}

zs = np.zeros(lmax_sinv + 1)
for k in cl_len_filt.keys():
    if k == "L":
        total_cl["L"] = cl_len_filt["L"]
    else:
        total_cl[k] = cl_len_filt.get(k, zs) + fgnd_cl.get(k, zs)

nl2d = pk.load(open(bdir + "inputs/nl2d/clnk150_dictfmt.pk", "rb"))
uK = core.G3Units.uK
nl2d = MapSpectraTEB(
    ffts=[
        nl2d["tfft"] * uK ** 2,
        nl2d["efft"] * uK ** 2,
        nl2d["bfft"] * uK ** 2,
    ],
    map_nx=npix,
    map_ny=npix,
    dx=reso_rad,
    dy=reso_rad,
)

sinvfilt150 = sl.cinv_utils.opfilt_teb.cl2sinv(
    total_cl, nl2d, tf150, nft=7.0, nfp=7.0, lmax=lmax_sinv
)

core.log_info("---- ivfs:cinv150", info="ParameterFile")
cinv150_len = sl.cinv.CinvFilt(
    obs150_len,
    sinvfilt150,
    ninvfilt150,
    bdir + "qest/%s/cinv150_len_t1p1/" % ivfs_prefix,
    eps_min=4.0e-4,
)
cinv150_len = sl.cinv.CinvFiltMasked(
    lxmin=lx_cut, lmax=lmax_cinv, cinv=cinv150_len
)

cinv150_unl = sl.cinv.CinvFilt(
    obs150_unl,
    sinvfilt150,
    ninvfilt150,
    bdir + "qest/%s/cinv150_unl_t1p1/" % ivfs_prefix,
    eps_min=4.0e-4,
)
cinv150_unl = sl.cinv.CinvFiltMasked(
    lxmin=lx_cut, lmax=lmax_cinv, cinv=cinv150_unl
)

cinv150_len_t2_nofg = sl.cinv.CinvFilt(
    obs150_len_t2_nofg,
    sinvfilt150,
    ninvfilt150,
    bdir + "qest/%s/cinv150_len_t2p1_nofg/" % ivfs_prefix,
    eps_min=4.0e-4,
)
cinv150_len_t2_nofg = sl.cinv.CinvFiltMasked(
    lxmin=lx_cut, lmax=lmax_cinv, cinv=cinv150_len_t2_nofg
)

cinv150_len_nofg = sl.cinv.CinvFilt(
    obs150_len_nofg,
    sinvfilt150,
    ninvfilt150,
    bdir + "qest/%s/cinv150_len_t1p1_nofg/" % ivfs_prefix,
    eps_min=4.0e-4,
)
cinv150_len_nofg = sl.cinv.CinvFiltMasked(
    lxmin=lx_cut, lmax=lmax_cinv, cinv=cinv150_len_nofg
)


# Separate these out, since different numbers of sims
# are used for different purposes

# everything, used for get_dat_teb()
ivflibs = [cinv150_len, cinv150_len_nofg, cinv150_len_t2_nofg]

ivflibs_mc_sims = [cinv150_len]  # evaluate for idxs in mc_sims
# evaluate for idxs in mc_sims_mf
ivflibs_mc_sims_n1 = [cinv150_len_nofg, cinv150_len_t2_nofg]
ivflibs_mc_sims_unl = [cinv150_unl]  # evaluate for idxs in mc_sims_unl


# Libraries we use to calculate quadratic estimates of phi, qest
core.log_info("---- qest", unit="ParameterFile")

# data (or sims treated identically to data)
qest_dd = sl.quadest.QuadEstLib(
    cl_unl,
    cl_len,
    cinv150_len,
    lib_dir=bdir + "qest/%s/par_%s/qest_len_dd/" % (ivfs_prefix, qest_prefix),
)

# data X sim.  Used for rdn0
qest_ds = sl.quadest.QuadEstLib(
    cl_unl,
    cl_len,
    cinv150_len,
    ivfs2=sl.cinv.CinvFiltData(cinv150_len),
    lib_dir=bdir + "qest/%s/par_%s/qest_len_ds/" % (ivfs_prefix, qest_prefix),
)

# simA X simB.  Used for N0
qest_ss = sl.quadest.QuadEstLib(
    cl_unl,
    cl_len,
    cinv150_len,
    ivfs2=sl.cinv.CinvFiltSim(cinv150_len, roll=2),
    lib_dir=bdir + "qest/%s/par_%s/qest_len_ss/" % (ivfs_prefix, qest_prefix),
)

# sim(t1p1) X sim(t2p1)
qest_ss2_nofg = sl.quadest.QuadEstLib(
    cl_unl,
    cl_len,
    cinv150_len_nofg,
    ivfs2=sl.cinv.CinvFiltSim(cinv150_len_t2_nofg, roll=0),
    lib_dir=bdir
    + "qest/%s/par_%s/qest_len_ss2_nofg/" % (ivfs_prefix, qest_prefix),
)

# simA X simB.  Used to subtract from N1
qest_ss_nofg = sl.quadest.QuadEstLib(
    cl_unl,
    cl_len,
    cinv150_len_nofg,
    ivfs2=sl.cinv.CinvFiltSim(cinv150_len_nofg, roll=2),
    lib_dir=bdir
    + "qest/%s/par_%s/qest_len_ss_nofg/" % (ivfs_prefix, qest_prefix),
)

# unlensed sims (unlensed input spectrum)
qest_uu = sl.quadest.QuadEstLib(
    cl_unl,
    cl_len,
    cinv150_unl,
    lib_dir=bdir + "qest/%s/par_%s/qest_unl_uu/" % (ivfs_prefix, qest_prefix),
)

qftkeys = ["ptt", "pee", "pte", "ptb", "peb", "p", "pp"]
qftlibs = [qest_dd, qest_ds, qest_ss, qest_uu, qest_ss_nofg, qest_ss2_nofg]


# qecl
core.log_info("---- qecl", unit="ParameterFile")
qecl_len_dd = sl.quadest_cl.QuadEstCl(
    qest_dd,
    lib_dir=bdir + "qest/%s/par_%s/qecl_len_dd/" % (ivfs_prefix, qest_prefix),
    mc_sims_mf=mc_sims_mf,
)
qecl_len_ds = sl.quadest_cl.QuadEstCl(
    qest_ds,
    lib_dir=bdir + "qest/%s/par_%s/qecl_len_ds/" % (ivfs_prefix, qest_prefix),
    mc_sims_mf=None,
)
qecl_len_ss = sl.quadest_cl.QuadEstCl(
    qest_ss,
    lib_dir=bdir + "qest/%s/par_%s/qecl_len_ss/" % (ivfs_prefix, qest_prefix),
    mc_sims_mf=None,
)
qecl_len_ss_nofg = sl.quadest_cl.QuadEstCl(
    qest_ss_nofg,
    lib_dir=bdir
    + "qest/%s/par_%s/qecl_len_ss_nofg/" % (ivfs_prefix, qest_prefix),
    mc_sims_mf=None,
)
qecl_len_ss2_nofg = sl.quadest_cl.QuadEstCl(
    qest_ss2_nofg,
    lib_dir=bdir
    + "qest/%s/par_%s/qecl_len_ss2_nofg/" % (ivfs_prefix, qest_prefix),
    mc_sims_mf=None,
)
qecl_len_uu = sl.quadest_cl.QuadEstCl(
    qest_uu,
    lib_dir=bdir + "qest/%s/par_%s/qecl_len_uu/" % (ivfs_prefix, qest_prefix),
    mc_sims_mf=mc_sims_mf,
)

qecl_len_dk = sl.quadest_cl.QuadEstCl(
    qest_dd,
    lib_dir=bdir + "qest/%s/par_%s/qxcl_len_dd/" % (ivfs_prefix, qest_prefix),
    qeB=sl.quadest.QuadEstLibKappa(qest_dd, scan.cmblib_proj_len_t1p1),
    mc_sims_mf=(mc_sims_mf, None),
)

qcllibs = [
    qecl_len_dd,
    qecl_len_ds,
    qecl_len_ss,
    qecl_len_uu,
    qecl_len_ss_nofg,
    qecl_len_ss2_nofg,
    qecl_len_dk,
]
qclkeys = list(zip(qftkeys, qftkeys)) + [
    ("ptt", "pte"),
    ("ptt", "pee"),
    ("ptt", "peb"),
    ("pee", "pte"),
    ("pte", "peb"),
]
