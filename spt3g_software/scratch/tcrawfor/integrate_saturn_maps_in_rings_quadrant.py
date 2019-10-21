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
p_main_220 = tctools.int_src(map220,4.,xcenter=180,ycenter=180,resolution=0.5)
p_main_90 = tctools.int_src(map90,4.,xcenter=180,ycenter=180,resolution=0.5)

# vector of radii 
rvec_inner = np.arange(25)*5. + 7.5
nrad = len(rvec_inner) - 1

npix_all_arr = np.zeros([5,nrad])
npix_used_arr = np.zeros([5,nrad])
pring150_arr = np.zeros([5,nrad])
pring220_arr = np.zeros([5,nrad])
pring90_arr = np.zeros([5,nrad])

masks = np.zeros([5,360,360])
for j in np.arange(4):
    q1 = j/2
    q2 = j - 2*q1
    masks[j,180*q1:180*(q1+1),180*q2:180*(q2+1)] = 1.
masks[4,:,:] = 1.

for i in np.arange(nrad):
    for j in np.arange(5):
        mask = masks[j,:,:]
        ptemp1 = tctools.int_src(map150_cut*mask,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
        ptemp2 = tctools.int_src(map150_cut*mask,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
        pring150_arr[j,i] = ptemp2 - ptemp1
        ptemp3 = tctools.int_src(map90_cut*mask,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
        ptemp4 = tctools.int_src(map90_cut*mask,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
        pring90_arr[j,i] = ptemp4 - ptemp3
        ptemp5 = tctools.int_src(cutmap*mask,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
        ptemp6 = tctools.int_src(cutmap*mask,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
        npix_used_arr[j,i] = ptemp6 - ptemp5
        ptemp7 = tctools.int_src(onesmap*mask,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
        ptemp8 = tctools.int_src(onesmap*mask,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
        npix_all_arr[j,i] = ptemp8 - ptemp7
        ptemp9 = tctools.int_src(map220_cut*mask,rvec_inner[i],xcenter=180,ycenter=180,resolution=0.5)
        ptemp10 = tctools.int_src(map220_cut*mask,rvec_inner[i+1],xcenter=180,ycenter=180,resolution=0.5)
        pring220_arr[j,i] = ptemp10 - ptemp9

pring90 = pring90_arr[4,:]*npix_all_arr[4,:]/npix_used_arr[4,:]
pring150 = pring150_arr[4,:]*npix_all_arr[4,:]/npix_used_arr[4,:]
pring220 = pring220_arr[4,:]*npix_all_arr[4,:]/npix_used_arr[4,:]
dpring90 = np.zeros(nrad)
dpring150 = np.zeros(nrad)
dpring220 = np.zeros(nrad)
for i in np.arange(nrad):
    temp1 = 4.*pring90_arr[0:4,i]*npix_all_arr[0:4,i]/npix_used_arr[0:4,i]
    pring90[i] = np.mean(temp1)
    dpring90[i] = np.std(temp1)/2.
    temp2 = 4.*pring150_arr[0:4,i]*npix_all_arr[0:4,i]/npix_used_arr[0:4,i]
    pring150[i] = np.mean(temp2)
    dpring150[i] = np.std(temp2)/2.
    temp3 = 4.*pring220_arr[0:4,i]*npix_all_arr[0:4,i]/npix_used_arr[0:4,i]
    pring220[i] = np.mean(temp3)
    dpring220[i] = np.std(temp3)/2.



