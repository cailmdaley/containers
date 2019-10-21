pro write_pointsources_sim09



ds = fltarr(6,3)

freqpairs = intarr(6,2)
freqpairs[0:2,0]=indgen(3)
freqpairs[0:2,1]=indgen(3)
freqpairs[3,0]=0
freqpairs[3,1]=1
freqpairs[4,0]=0
freqpairs[4,1]=2
freqpairs[5,0]=1
freqpairs[5,1]=2


;radio
fr=[95,151.7,219.3]
alphar = -0.6
dr=1.3
for i=0,5 do ds[i,0] = sqrt($
                             scale_ps(dr,fr[1],fr[freqpairs[i,0]],alphar) *$
                             scale_ps(dr,fr[1],fr[freqpairs[i,1]],alphar))

;dsfg
f=[95,154.2,221.4]
alpha=3.6
dp = 7.7
for i=0,5 do ds[i,1] = sqrt($
                             scale_ps(dp,f[1],f[freqpairs[i,0]],alpha) *$
                             scale_ps(dp,f[1],f[freqpairs[i,1]],alpha))

;note this is ell-like
dc = 5.9
for i=0,5 do ds[i,2] = sqrt($
                             scale_ps(dc,f[1],f[freqpairs[i,0]],alpha) *$
                             scale_ps(dc,f[1],f[freqpairs[i,1]],alpha))

comb = total(ds[*,0:1],2)

freqs=['90','150','220']

ofile='/home/cr/cmb_models/foreground_sim09_'+freqs+'.txt'
l=lindgen(25000)+2
psfac = l*(l+1d)/(2*!pi) / ((3000.*3001.)/(2.*!pi))
clusfac = l/3000.

for i=0,2 do begin
    openw,1,ofile[i]
    
    dl = psfac*comb[i] + clusfac*ds[i,2]
    for j=0,24999 do $
      printf,1,l[j],dl[j]
    
    
    close,1

endfor


end
