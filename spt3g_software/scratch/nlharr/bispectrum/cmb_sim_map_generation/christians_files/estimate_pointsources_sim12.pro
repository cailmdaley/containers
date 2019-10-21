function estimate_pointsources_sim12



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
fr=[97.3,150.2,214.1]
alphar = -0.53
dr=1.28
for i=0,5 do ds[i,0] = sqrt($
                             scale_ps(dr,fr[1],fr[freqpairs[i,0]],alphar) *$
                             scale_ps(dr,fr[1],fr[freqpairs[i,1]],alphar))

;dsfg
f=[97.9, 153.8, 219.6]
alpha=3.56
dp = 7.54
for i=0,5 do ds[i,1] = sqrt($
                             scale_ps(dp,f[1],f[freqpairs[i,0]],alpha) *$
                             scale_ps(dp,f[1],f[freqpairs[i,1]],alpha))

;note this is ell^0.8-like
dc = 6.25
for i=0,5 do ds[i,2] = sqrt($
                             scale_ps(dc,f[1],f[freqpairs[i,0]],alpha) *$
                             scale_ps(dc,f[1],f[freqpairs[i,1]],alpha))

comb = total(ds[*,0:1],2)
return,{ps:comb,dsfg:reform(ds[*,1]),radio:reform(ds[*,0]),$
        ellclus:reform(ds[*,2])}
end
