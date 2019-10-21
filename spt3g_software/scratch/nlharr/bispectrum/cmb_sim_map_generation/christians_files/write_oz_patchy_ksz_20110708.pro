pro write_oz_patchy_ksz_20110708

a = read_ascii('~cr/cmb_models/patchy_ksz_oz_cl768_750Mpc_eff20_deff0.cl')


inl=[-1,2,reform(a.field1[0,*])]
incl = [0,0,reform(a.field1[1,*])]

ll=dindgen(11001)
tdl = interpol(incl,inl,ll)
dl = smooth(tdl,51) * 1e12 * 2.726^2
cl = dl /( ll * (ll+1.)/2./!dpi)
cl[0:2]=0
;also make one without the lowell
dl2=dl
dlb=dl2[500]
dl2[0:500]=dindgen(501)*dl[500]/500.
cl2 = dl2 /( ll * (ll+1.)/2./!dpi)
cl2[0:2]=0



plot,ll,dl,xr=[0,10000]
oplot,ll,dl2,color=200

clfile='~cr/paramfits/cosmomc.r11/ptsrc/cl_ksz_oz_patchy_20110708.txt'
dlfile='~cr/paramfits/cosmomc.r11/ptsrc/dl_ksz_oz_patchy_20110708.txt'
cl2file='~cr/paramfits/cosmomc.r11/ptsrc/cl_ksz_oz_patchy_nolowell_20110708.txt'
dl2file='~cr/paramfits/cosmomc.r11/ptsrc/dl_ksz_oz_patchy_nolowell_20110708.txt'
openw,1,clfile
openw,2,dlfile
for i=0,10000 do begin
    printf,1,format='(i,e)',ll[i],cl[i]
    printf,2,format='(i,e)',ll[i],dl[i]
endfor
close,1
close,2
openw,1,cl2file
openw,2,dl2file
for i=0,10000 do begin
    printf,1,format='(i,e)',ll[i],cl2[i]
    printf,2,format='(i,e)',ll[i],dl2[i]
endfor
close,1
close,2
stop

end
