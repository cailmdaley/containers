pro sim_cluster_distort


n1=7040*2
n2=3944*2
res=0.125

xs=[3620,1400,700,3420,5550,6200]*2
ys=[1872,3350,1400,2072,3350,1350]*2
thetas=[6.,6,6,2,2,2.]

nclus = 6

maxr = 10*max(thetas)+3
nr = 2*maxr + 1
xx = (findgen(nr)-maxr) # (fltarr(nr)+1)
yy = transpose(xx)
distsq = xx*xx + yy*yy

map = fltarr(n1,n2)

for i=0,nclus -1 do begin
    
    dmax = 10*thetas[i]
    rsq = thetas[i]^2
    
    inds = where(distsq gt dmax*dmax)
    tmap = -1./(1 + distsq/rsq)
    tmap(inds) = 0.0
    tt = total(tmap)
    tmap *= 100./tt
    
    map[xs[i]-maxr:xs[i]+maxr,$
        ys[i]-maxr:ys[i]+maxr]=tmap

endfor

ell = findgen(30001)
bl = calc_bl(ell,1.,beamfile='/home/cr/data/spt_beams_2009/bl_spt_2009_150.txt')
qvec=ell/(60.*180./!pi)

omap = apply_1d_filter_to_map(map,.125,.125,qvec,bl)


openw,1,'/data/cr/lotsomaps_clustest.dat'
writeu,1,float(omap)
close,1

end
