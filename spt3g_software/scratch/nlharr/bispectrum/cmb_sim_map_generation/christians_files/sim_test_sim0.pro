function beam_l,sigma_b,l
a=0.85
return,  (1.-a) * exp(-.5 * (.00292*l)^1.8) + a * exp(-.5 * (sigma_b*l)^1.5)
end
function beam_adj,dsigmas,fri,frj,l
 sigma_b0 = [0.00017172336, 0.00016751516]
 sigma_b = sigma_b0 + dsigmas
 
 return,  beam_l(sigma_b0(fri),l)*beam_l(sigma_b0(frj),l)/beam_l(sigma_b(fri),l)/beam_l(sigma_b(frj),l)
 
end

pro sim_mf_diff_s10

nsim = 1

seed=7
;

;    ! Effective frequencies for the bands                                                             
;   ! The order is (dutsy clustered,
                                ;   dusty poisson, radio poisson, ksz,
                                ;   tsz)                          
;  real, dimension(5,2), parameter :: spt_eff_fr \
   ;     = (/ (/ 154.2, 154.2, 151.7, 150, 153/), \
 ;            (/ 221.4, 221.4, 219.3, 220, 220.2/) /)


l = dindgen(9951)+50


readcol,'/home/cr/paramfits/cosmomc.s10/camb/ml_20100826/ml_l10000_acc2_lensedCls.dat',l2,dl2,t1,t2,t3
dlcmb = [[dl2[l-2]],[dl2[l-2]],[dl2[l-2]]]

llfac0 = l*(l+1d)/(2d*!dpi)
llfac = llfac0 / ((3000* 3001d)/(2d*!dpi))
dlradio150 = 1.28*llfac
dlradio220 = 0.95*llfac
dlradio =[[dlradio150],[sqrt(dlradio150*dlradio220)],[dlradio220]]

dldsfg = [[llfac0*6.6e-6],[llfac0* 2.24439e-05],[llfac0*7.73299e-05 ]]

readcol,'/home/cr/paramfits/cosmomc.s10/ptsrc/tsz_sehgal.dat',l2,dtsz
dtszs = dlcmb*0.0
dtszs[*,0]=[fltarr(50),dtsz*llfac0]*0.5

readcol,'/home/cr/paramfits/cosmomc.s10/ptsrc/ksz_sehgal.dat',l2,dksz
dkszs = dlcmb*0.0
dkszs[*,0]=[dksz[48:*]*llfac0]
dkszs[*,1]=[dksz[48:*]*llfac0]
dkszs[*,2]=[dksz[48:*]*llfac0]

readcol,'/home/cr/paramfits/cosmomc.s10/ptsrc/clustered_150.dat',l2,cclus
dclus = dlcmb*0.0
dlclus = cclus[48:*]*llfac0
dlclus /=dlclus[2950]
dlclus *=4.5
temp220 = (l/3000.)^(-.0000184825*2.*(221.4-154.2))
dlclus220 = dlclus * temp220 * (221.4/154.2)^(2*3.88)
dclus[*,0]=dlclus
dclus[*,1]=sqrt(dlclus220*dlclus)
dclus[*,2]=dlclus220

dlth = dlradio + dkszs + dtszs + dlcmb + dldsfg + dclus
;stop

wins = fltarr(9951,45)
for i=1,45 do begin
    readcol,'~cr/s10/window_'+strtrim(string(i),2),l3,w3
    wins[*,i-1]=w3[0:9950]
endfor
dltmp = fltarr(9951,3)

readcol,'~cr/s10/covariance_150_cross_220_s10.txt',cov
cov = reform(cov,45,45)
evec=cov
TRIRED, evec, D, E  
; Compute the eigenvalues (returned in vector D) and  
; the eigenvectors (returned in the rows of the array A):  
TRIQL, D, E, evec  
err = sqrt(d)

banddef=[1000, 2000, 2200, 2400, 2600, 2800, 3000, 3400, 3800, 4200, 4600, $
            5000, 5900, 6800, 7700, 8600, 9500]
center=([0, banddef]+banddef)/2
i0=2
i1=16
banddef=banddef[i0:i1]
center=center[i0:i1]
nii=i1-i0+1

alpha=-0.325
arunname0=string(alpha)

caluncert150=0.036
 ; 3.6% + 5% systematic uncertainty on calibration
caluncertX=sqrt(0.036^2+0.05^2)
mapcaluncert=sqrt(caluncert150^2+(alpha/(alpha+1.0))^2*caluncertX^2)

speccaluncert=(1+mapcaluncert)^2-1

berr=fit_gauss_err(alpha=-1*alpha, /nongauss,ell_in=center)


chisq = fltarr(nsim)
for i=0,nsim-1 do begin
    a150 = 1 + randomn(seed) *.036
    rX = 1 + randomn(seed)*.062
    
    calfacs = a150^2 * [1.,rx,rx*rx]
    
    dsigma150 = randomn(seed)*2.86e-6
    dsigma220 = randomn(seed)*6.83e-6
    
    dsigmas = [dsigma150,dsigma220]
  

    dltmp[*,0] = dlth[*,0] * calfacs[0]* beam_adj(dsigmas,0,0,l)
    dltmp[*,1] = dlth[*,1] * calfacs[1]* beam_adj(dsigmas,0,1,l)
    dltmp[*,2] = dlth[*,2] * calfacs[2]* beam_adj(dsigmas,1,1,l)
 
    print,i,calfacs,dsigmas
    ;now put into window functions!
    continue
    
    bs = transpose(wins[*,0:14]) # dltmp
    bs[*,1] = transpose(wins[*,15:29]) # dltmp[*,1]
    bs[*,2] = transpose(wins[*,30:44]) # dltmp[*,2]
    bs = reform(bs,45)
    
    ;add errors
    berrdiag =  randomn(seed,45)*err
    berr = evec * berrdiag
    
    chisq[i] = berr # invert(cov) # berr

    bsfinal = bs + berr

;    continue
;    print,'shouldnt be here'
;    stop

    ; now write to disk
    openw,lun,/get_lun,'~cr/s10_mf_test/test_'+strtrim(string(i),2)
    for j=0,44 do $
      printf,lun,j mod 15,bsfinal[j]
    free_lun,lun
    

    
    ;now do MF -> differenced spectrum
    release=spectrum_subtraction(banddef, bsfinal, cov, $
                                 wins, alpha=alpha)
    release.bandpowers/=(1+alpha)^2 * 1e12
    release.covariance/=(1+alpha)^4 * 1e24
    release.windowfunctions/=(1+alpha)^2
    
    cmcfile='~cr/s10_mf_test/alpha0p325/Spectrum_2008_alpha0p325_'+strtrim(string(i),2)+'.newdat'
    
    

    if i eq 0 then $
      spectrum_to_oz_cosmomc, banddef[*], $
      release.bandpowers[*], $
      release.covariance[*, *], cmcfile=cmcfile, $
      window=release.windowfunctions, $
      winfuncpref='window_sim/window_', $
      comment=comment, caluncert=speccaluncert, $
      fullbeamerr=berr[*], ellmin=2000 $
    else $
      spectrum_to_oz_cosmomc, banddef[*], $
      release.bandpowers[*], $
      release.covariance[*, *], cmcfile=cmcfile, $
      winfuncpref='window_sim/window_', $
      comment=comment, caluncert=speccaluncert, $
      fullbeamerr=berr[*], ellmin=2000 
    
    
endfor

print,moment(chisq)
end
