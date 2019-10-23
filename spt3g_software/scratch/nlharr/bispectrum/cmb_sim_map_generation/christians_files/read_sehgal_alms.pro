function read_sehgal_alms,file
;this pro fails to work!
stop


a=READfits(file,EXTEN_NO=1)
index = long64(256)*4*a[0,*]+$
  long64(256)*2*a[1,*]+$
  long64(256)*a[2,*]+$
  a[3,*]-1
;index = l*(l+1) + m + 1 - 1
openw,1,'tmp'
writeu,1,reverse(a[4:7,*],1)
writeu,1,reverse(a[8:11,*],1)
close,1
nn=n_elements(index)
ireal=fltarr(nn)
iimag=fltarr(nn)
openr,1,'tmp'
readu,1,ireal
readu,1,iimag
close,1
stop
openw,1,'tmp'
writeu,1,reverse(a[4:7,0])
writeu,1,reverse(a[4:7,1])
writeu,1,reverse(a[4:7,2])
close,1
openr,1,'tmp'
b=fltarr(3)
readu,1,b
close,1
stop
alms=ireal+complex(0,1)*iimag
stop
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
