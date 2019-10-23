pro convert_sehgal_to_text,file,ofile
;takes a anafast spectrum based on a Sehgal healpix map

;file='/home/cr/data/sehgal_model/spectra_148_ksz.fits'
;fxbopen,lun,file,1
;        fxbreadm,lun,[1],cl
;        fxbclose,lun

;fxopen,lun,file,1
;        fxread,file,cl,header,extension=1
;        fxclose,lun

a=READfits(file,EXTEN_NO=1)

nl=n_elements(a[0,*])
l = findgen(nl)

cl = dblarr(nl)
for i=0,nl-1 do begin
    cl[i]=double(string(a[*,i]))

endfor
conversion = (2.726e6/1.07248e9)^2 ; from Jy2/sr to muK2
idl = cl*l*(l+1)/(2.*!pi)*conversion
;*(2.73e-3)^2
;added a smoothing to reduce noise)
dl = smooth(idl,31)

;ofile='~cr/paramfits/cosmobar/data/dl_ksz_sehgal.txt'
openw,lun,ofile,/get_lun
for i=0,nl-1 do $
  printf,lun,format='(i,e)',l[i],dl[i]

free_lun,lun
;stop

end
