pro write_oz_patchy_ksz_20110503

a = read_ascii('~cr/cmb_models/patchy_kSZ_res768.cl')
inl=[-1,2,reform(a.field1[0,*])]
incl = [0,0,reform(a.field1[1,*])]

ll=dindgen(11001)
tdl = interpol(incl,inl,ll)
dl = smooth(tdl,51) * 1e12 * 2.726^2
cl = dl /( ll * (ll+1.)/2./!dpi)
cl[0:2]=0


plot,ll,dl,xr=[0,10000]

clfile='~cr/paramfits/cosmomc.r11/ptsrc/cl_ksz_oz_patchy_20110503.txt'
dlfile='~cr/paramfits/cosmomc.r11/ptsrc/dl_ksz_oz_patchy_20110503.txt'
openw,1,clfile
openw,2,dlfile
for i=0,10000 do begin
    printf,1,format='(i,e)',ll[i],cl[i]
    printf,2,format='(i,e)',ll[i],dl[i]
endfor
close,1
close,2
stop

end
