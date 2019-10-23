function estimate_pointsources_sim12_spire

inmjy = [0,0,0,1,1,1]
f=[97.9, 153.8, 219.6,600,857.,1200.]
fr=[97.3,150.2,214.1,600,857.,1200.]

ds = fltarr(21,3)

freqpairs = intarr(21,2)
k=0
for i=0,5 do for j=i,5 do begin
    freqpairs[k,0]=i
    freqpairs[k,1]=j
k++
endfor
;old ordering
;freqpairs[0:2,0]=indgen(3)
;freqpairs[0:2,1]=indgen(3)
;freqpairs[3,0]=0
;freqpairs[3,1]=1
;freqpairs[4,0]=0
;freqpairs[4,1]=2
;freqpairs[5,0]=1
;freqpairs[5,1]=2


;radio

alphar = -0.53
sigmarsq = 0.1
dr=1.28
for i=0,20 do ds[i,0] = $
  scale_ps_wscatter(dr,fr[1],fr[freqpairs[i,0]],fr[freqpairs[i,1]],alphar,sigmarsq)


;dsfg

t=12.
beta=2.
sigmasq=0.1
dp = 7.54
for i=0,20 do ds[i,1] = $
  scale_ps_modbb(dp,fr[1],fr[freqpairs[i,0]],fr[freqpairs[i,1]],t,beta,sigmasq)


;note this is ell^0.8-like
dc = 6.25
t=12.
beta=2.
for i=0,20 do ds[i,2] = $
  scale_ps_modbb(dc,fr[1],fr[freqpairs[i,0]],fr[freqpairs[i,1]],t,beta,sigmasq)



for i=0,20 do begin
    if inmjy[freqpairs[i,0]] then begin
        ds[i,1:2] *= mjyarcmin_to_ukarcmin(1.,freq=f[freqpairs[i,0]],/invert)
        ds[i,0] *= mjyarcmin_to_ukarcmin(1.,freq=fr[freqpairs[i,0]],/invert)
    endif
    if inmjy[freqpairs[i,1]] then begin
        ds[i,1:2] *= mjyarcmin_to_ukarcmin(1.,freq=f[freqpairs[i,1]],/invert)
        ds[i,0] *= mjyarcmin_to_ukarcmin(1.,freq=fr[freqpairs[i,1]],/invert)
    endif
endfor

comb = total(ds[*,0:1],2)
return,{ps:comb,dsfg:reform(ds[*,1]),radio:reform(ds[*,0]),$
        ellclus:reform(ds[*,2]),fr:fr,fdust:f,pairs:freqpairs}
end
