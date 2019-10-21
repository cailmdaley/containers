# assume map150 and map90 exist and are 360x360 with 0.5' resolution

# get rid of poly trench
map150_cut=map150.copy()
map150_cut[160:200,:]=0.
map90_cut=map90.copy()
map90_cut[160:200,:]=0.
map220_cut=map220.copy()
map220_cut[160:200,:]=0.
cutmap=np.zeros([360,360])+1.
cutmap[160:200,:]=0.
onesmap=np.zeros([360,360])+1.

# main lobe power
p_main_150 = tctools.int_src(map150,4.,xcenter=180,ycenter=180,resolution=0.5)
p_main_90 = tctools.int_src(map90,4.,xcenter=180,ycenter=180,resolution=0.5)
p_main_220 = tctools.int_src(map220,4.,xcenter=180,ycenter=180,resolution=0.5)

# vector of radii 
rvec_inner = np.arange(25)*5. + 7.5
nrad = len(rvec_inner) - 1

npix_all = np.zeros(nrad)
npix_used = np.zeros(nrad)
pring150_orig = np.zeros(nrad)
pring90_orig = np.zeros(nrad)
pring220_orig = np.zeros(nrad)

for i in np.arange(nrad):
    ptemp1 = tctools.int_src(map150_cut,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
    ptemp2 = tctools.int_src(map150_cut,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
    pring150_orig[i] = ptemp2 - ptemp1
    ptemp3 = tctools.int_src(map90_cut,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
    ptemp4 = tctools.int_src(map90_cut,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
    pring90_orig[i] = ptemp4 - ptemp3
    ptemp5 = tctools.int_src(cutmap,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
    ptemp6 = tctools.int_src(cutmap,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
    npix_used[i] = ptemp6 - ptemp5
    ptemp7 = tctools.int_src(onesmap,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
    ptemp8 = tctools.int_src(onesmap,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
    npix_all[i] = ptemp8 - ptemp7
    ptemp9 = tctools.int_src(map220_cut,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
    ptemp10 = tctools.int_src(map220_cut,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
    pring220_orig[i] = ptemp10 - ptemp9
