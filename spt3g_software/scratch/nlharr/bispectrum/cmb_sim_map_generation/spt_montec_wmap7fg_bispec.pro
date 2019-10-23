pro spt_montec_wmap7fg_bispec, freq, n1, n2, n, scmb, sps, sps2, stsz, output_file, reso_arcmin, nrealizations
;n1=3228
;n2=2976
;n=3240 ; the square size of the map
;scmb=105563
;sps=99109
;sps2=32457277
;stsz = 15486277
;output_file = '/data30/nlharr/Bispectrum2500D/sims/test_file'
;reso_arcmin = 0.5
;nrealizations = 100




;; In python maps = np.fromfile(output_file, dtype='float32').reshape((nsim,n2,n1))

pi=acos(-1.d0)
resolution=(reso_arcmin/60d0)/180d0*pi

modelling_folder = "/data30/nlharr/Bispectrum2500D/modelling_files/"

;element == 0 cmb
;element == 1 SZ
;element == 2 PS

;waka=20
waka=pi/(2d0*n*resolution)


;pi/waka=size of the map
;only l>2*waka part of the fluctuation are present
;map resolution= pi/(waka*2*n), 

df=float(waka)/pi

dtheta=1./df/2./n
;frequency resolution
;frequency span= -df*n ~ df*n (our map is 2n by 2n)

spectrum05=complexarr(n*2,n*2) ; will be combined for all years

gar=intarr(2)


calc=dblarr(25000)
calccmb=dblarr(25000)
ll=dindgen(25000)+2
dlfac = ll*(ll+1.)/(2*!pi)
fwhm2sig = 1./sqrt(8.*alog(2.))

;2009 used these seeds:   
;scmb=11 
;sps=23
;sps2=31
;2008 is using these:

;compstr='ps'
;calc = spt_th_cl_run2(ll,150,/onlyps)
calc = dlfac / ((3000.*3001.)/(2.*!pi))
                                ; I've normalized with respect to
                                ; D3000=1
calcell = (ll/3000.)^0.8


facs = estimate_pointsources_sim12_spire()
matrix = dblarr(6,6)
ellmatrix=matrix
k=0
for i=0,5 do for j=i,5 do begin
    matrix[i,j]=facs.ps[k]
    matrix[j,i]=facs.ps[k]
    ellmatrix[i,j]=facs.ellclus[k]
    ellmatrix[j,i]=facs.ellclus[k]
    k++
endfor
evecs=matrix
trired,evecs,evals,e
triql,evals,e,evecs

ii0=1
if fix(freq) eq 220 then ii0=2
if fix(freq) eq 90 then ii0=0
if fix(freq) eq 600 then ii0=3
if fix(freq) eq 857 then ii0=4
if fix(freq) eq 1200 then ii0=5

amps = sqrt(evals)
effs=fltarr(6)
effs[0] = amps[0]*evecs[ii0,0]
effs[1] = amps[1]*evecs[ii0,1]
effs[2] = amps[2]*evecs[ii0,2]
effs[3] = amps[3]*evecs[ii0,3]
effs[4] = amps[4]*evecs[ii0,4]
effs[5] = amps[5]*evecs[ii0,5]

;check if planned N makes sense
jj=abs(evals)
ind = where(jj gt 1e-4*max(jj),nbig,complement=cnd)
if nbig ne 4 then stop
if nbig lt 6 then effs[cnd]=0

evecs=ellmatrix
trired,evecs,evals,e
triql,evals,e,evecs

jj=abs(evals)
ind = where(jj gt 1e-4*max(jj),nbig,complement=cnd)
if nbig ne 3 then stop
ampsclus = sqrt(evals)

effsclus=fltarr(6)
effsclus[0] = ampsclus[0]*evecs[ii0,0]
effsclus[1] = ampsclus[1]*evecs[ii0,1]
effsclus[2] = ampsclus[2]*evecs[ii0,2]
effsclus[3] = ampsclus[3]*evecs[ii0,3]
effsclus[4] = ampsclus[4]*evecs[ii0,4]
effsclus[5] = ampsclus[5]*evecs[ii0,5]

if nbig lt 6 then effsclus[cnd]=0

tcmb=float(2.726)
if ii0 le 2 then begin
    calccmb= spt_th_cl_run3_2009(ll, freq, /onlycmb)
    
;also read in kSZ
    a = read_ascii(modelling_folder + 'dl_ksz_CSFplusPATCHY_13sep2011_norm1_fake25000.txt')
    ksz=ll*0.0
    loff = a.field1[0,0]
    ksz[0:24998]=a.field1[1,(2-loff):24998+(2-loff)]
;2.0 uk2 at ell=3000
    ksz *= 2.0
    
    calccmb +=ksz
    
;tSZ version
    a = read_ascii(modelling_folder + 'dl_shaw_tsz_s10_153ghz_fake25000.txt')
    tsz = ll*0.0
    loff = a.field1[0,0]
    tsz[0:24998]=a.field1[1,(2-loff):24998+(2-loff)]
    fszs=[97.6,152.9,218.1]
    scale = (fxsz(fszs[ii0])/fxsz(153.))^2
    tsz *= scale
;stop
    calccmb=(calccmb/1d12)/tcmb^2
    calctsz = (tsz/1d12)/tcmb^2
endif

calc=(calc/1d12)/tcmb^2
calcell=(calcell/1d12)/tcmb^2


cmb05=fltarr(25001)
cmb05b=fltarr(25001)
cmb05c=fltarr(25001)
cmb05d=fltarr(25001)
;s=8
;s=7
;s=6


beamsize05=1.19
if keyword_set(nofilt) then beamsize05 = sqrt(0.03*0.03)
fwhm05=beamsize05*pi/60./180.
sig05=fwhm05/sqrt(8.*alog(2.))
sigl05=1./sig05
;this part just for RMS calculation. The power spectrum
;has not actually been convolved

fstr = strtrim(string(fix(freq)),2)


; just setting the beam to be 1
bl = ll * 0.0 + 1.0

;if ii0 le 2 then begin
;    bdir = '/home/cr/data/spt_beams_2010/v1/'
;    fbeam = 'blgrid_2010_'+fstr+'.txt'
;    bl = calc_bl(ll,1.0,beamfile=bdir+fbeam)
;endif else begin
;    fwhm=[35.,24.,18.]/60.
;    bl = calc_bl(ll,fwhm[ii0-3])
;endelse

pix = resolution;

var1=0.
var2=0.
for l=2, 25000 do begin
    cmb05[l]=(2.*pi*calc[l-2])/(l*(l+1.)) * bl[l-2]^2 ;*exp(-float(l)^2/sigl05^2) / $
    cmb05c[l]=(2.*pi*calcell[l-2])/(l*(l+1.)) * bl[l-2]^2 ;*exp(-float(l)^2/sigl05^2) / $

    var1=var1+(2*l+1.)*cmb05[l]/4./pi*exp(-float(l)^2/sigl05^2)

endfor
if ii0 le 2 then begin
    for l=2, 25000 do begin
        cmb05b[l]=(2.*pi*calccmb[l-2])/(l*(l+1.)) * bl[l-2]^2
        cmb05d[l]=(2.*pi*calctsz[l-2])/(l*(l+1.)) * bl[l-2]^2
        var2+=(2*l+1.)*cmb05b[l]/4./pi*exp(-float(l)^2/sigl05^2)
    endfor
    print,'CMB Standard deviation =',2.726e3*sqrt(var2),'  mK'
    cmb05b[0:1]=0.
    cmb05d[0:1]=0.
endif 

cmb05c[0:100]=cmb05c[101]
;stop
print,'PS Standard deviation =',2.726e3*sqrt(var1),'  mK'

;stop
var2=0.

cmb05[0:1]=0.
cmb05c[0:1]=0.


print, output_file
;stop
openw,5,output_file

ddd=dist(2*n,2*n)*2.*waka
lll=ddd[0:n,*]
lllinds = where(lll gt 0 and lll lt 25000,nkeep)
amps = dblarr(n+1,2*n)
amps(lllinds) = sqrt((df^2.*cmb05[fix(lll(lllinds))])/2.)
ampsc = dblarr(n+1,2*n)
ampsc(lllinds) = sqrt((df^2.*cmb05c[fix(lll(lllinds))])/2.)

if ii0 le 2 then begin
    ampsb = dblarr(n+1,2*n)
    ampsb(lllinds) = sqrt((df^2.*cmb05b[fix(lll(lllinds))])/2.)
    ampsd = dblarr(n+1,2*n)
    ampsd(lllinds) = sqrt((df^2.*cmb05d[fix(lll(lllinds))])/2.)
endif

rtmp1=complexarr(n+1,2*n)
rtmp2=complexarr(n+1,2*n)
rtmp3=complexarr(n+1,2*n)
rtmp4=complexarr(n+1,2*n)

for realization=0,nrealizations-1 do begin
;for realization=0,1-1 do begin
    print,realization
    
    rtmp1(lllinds)=amps(lllinds)*complex(randomn(sps,nkeep),randomn(sps,nkeep))
    rtmp2(lllinds)=amps(lllinds)*complex(randomn(sps,nkeep),randomn(sps,nkeep))
    rtmp3(lllinds)=amps(lllinds)*complex(randomn(sps,nkeep),randomn(sps,nkeep))
    rtmp4(lllinds)=amps(lllinds)*complex(randomn(sps,nkeep),randomn(sps,nkeep))
    spectrum05(0:n,0:2*n-1)=rtmp1*effs[2]+rtmp2*effs[3]+rtmp3*effs[4]+rtmp4*effs[5]


    ;now add ell-clustered
    rtmp1(lllinds)=ampsc(lllinds)*complex(randomn(sps2,nkeep),randomn(sps2,nkeep))
    rtmp2(lllinds)=ampsc(lllinds)*complex(randomn(sps2,nkeep),randomn(sps2,nkeep))
    rtmp3(lllinds)=ampsc(lllinds)*complex(randomn(sps2,nkeep),randomn(sps2,nkeep))
    spectrum05(0:n,0:2*n-1)+=rtmp1*effsclus[3]+rtmp2*effsclus[4]+rtmp3*effsclus[5]

    ;now add CMB & tsz
    if ii0 le 2 then begin
        ;CMB
        rtmp1(lllinds)=ampsb(lllinds)*complex(randomn(scmb,nkeep),randomn(scmb,nkeep))
        spectrum05(0:n,0:2*n-1)+=rtmp1
                                ;now add tSZ
        rtmp3(lllinds)=ampsd(lllinds)*complex(randomn(stsz,nkeep),randomn(stsz,nkeep))
        spectrum05(0:n,0:2*n-1)+=rtmp3
    endif
    

;    if finite(max(spectrum05)) eq 0 then stop
;Gaussian statistics
;
;	e s s s o o o o
;	e s s s o o o o
;	e s s s o o o o
;	y s s s x o o o
;	s s s s s o o o
;	s s s s s o o o
;	s s s s s o o o
;	r s s s z q q q
;
;	s: specified complex points
;	x,y,z,r: specified real points
;	o: calculated from symmetry about x (in order for map to be real)
;	e: calculated from symmetry about y
;	q: calculated from symmetry about z

; now for the "o" part:

jj=1+indgen(2*n-1)
for i=1, n-1 do begin
    spectrum05(2*n-i,2*n-jj) = conj(spectrum05(i,jj))
endfor

jj = n+1+indgen(n-1)
spectrum05(n,jj)=conj(spectrum05(n,2*n-jj))
spectrum05(0,jj)=conj(spectrum05(0,2*n-jj))
spectrum05(jj,0)=conj(spectrum05(2*n-jj,0))

; a vertical "o" string is missed:
;spectrum(n,j)=conj(spectrum(n,2*n-j))
; "e":
;spectrum(0,j)=conj(spectrum(0,2*n-j))
; "q":
;spectrum(j,0)=conj(spectrum(2*n-j,0))


spectrum05[0,0]=complex(spectrum05[0,0],0)
spectrum05[n,n]=complex(spectrum05[n,n],0)
spectrum05[n,0]=complex(spectrum05[n,0],0)
spectrum05[0,n]=complex(spectrum05[0,n],0)

map05=tcmb*float(fft(spectrum05,1))

;stop
writeu,5,float(map05[0:n1-1,0:n2-1])

endfor

close,5
close,/all

;stop

end




