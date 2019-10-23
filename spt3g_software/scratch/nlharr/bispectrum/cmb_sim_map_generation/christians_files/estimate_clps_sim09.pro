function estimate_clps_sim09

uktojy = [2.03e8,3.99e8,4.84e8]
freqs=[90.,150.,220.]

;clps, alpha90-150; alpha 150-220; sigma_alpha^2
radiopars = [9.1e-7,-0.1,-0.5,0.04]
;radiopars = [9.1e-7,-0.5,-0.5,0.04]
dsfgpars  = [6.6e-6, 3.7,3.7,0.09]

clpss = fltarr(6,2) 
freqpairs = clpss
freqpairs[0:2,0]=freqs
freqpairs[0:2,1]=freqs
freqpairs[3,0]=90.
freqpairs[3,1]=150.
freqpairs[4,0]=90.
freqpairs[4,1]=220.
freqpairs[5,0]=150.
freqpairs[5,1]=220.

pars = radiopars
clpss[0,0] = pars[0] * (uktojy[1]^2 / uktojy[0]^2 ) * $
  (freqs[0]/freqs[1])^(2.*pars[1]) * ( 1. + pars[3]*(alog(freqs[0]/freqs[1]))^2)
clpss[1,0] = pars[0] 
clpss[2,0] = pars[0] * (uktojy[1]^2 / uktojy[2]^2 ) * $
  (freqs[2]/freqs[1])^(2.*pars[2]) * ( 1. + pars[3]*(alog(freqs[2]/freqs[1]))^2)
clpss[3,0] = pars[0] * (uktojy[1] / uktojy[0] ) * $
  (freqs[0]/freqs[1])^(pars[1]) 
clpss[4,0] = pars[0] * (uktojy[1]^2 / (uktojy[0]*uktojy[2]) ) * $
  (freqs[0]/freqs[1])^(pars[1]) * (FREQS[2]/freqs[1])^(pars[2]) 
clpss[5,0] = pars[0] * (uktojy[1] / uktojy[2] ) * $
  (freqs[2]/freqs[1])^(pars[2]) 

pars = DSFGpars
clpss[0,1] = pars[0] * (uktojy[1]^2 / uktojy[0]^2 ) * $
  (freqs[0]/freqs[1])^(2.*pars[1]) * ( 1. + pars[3]*(alog(freqs[0]/freqs[1]))^2)
clpss[1,1] = pars[0] 
clpss[2,1] = pars[0] * (uktojy[1]^2 / uktojy[2]^2 ) * $
  (freqs[2]/freqs[1])^(2.*pars[2]) * ( 1. + pars[3]*(alog(freqs[2]/freqs[1]))^2)
clpss[3,1] = pars[0] * (uktojy[1] / uktojy[0] ) * $
  (freqs[0]/freqs[1])^(pars[1]) 
clpss[4,1] = pars[0] * (uktojy[1]^2 / (uktojy[0]*uktojy[2]) ) * $
  (freqs[0]/freqs[1])^pars[1] * (FREQS[2]/freqs[1])^(pars[2]) 
clpss[5,1] = pars[0] * (uktojy[1] / uktojy[2] ) * $
  (freqs[2]/freqs[1])^(pars[2]) 


comb = total(clpss,2)
return,{cls:comb,dsfg:reform(clpss[*,1]),radio:reform(clpss[*,0])}

end
