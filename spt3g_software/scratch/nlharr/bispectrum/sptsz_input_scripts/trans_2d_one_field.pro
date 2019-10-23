
;example file
;trans_2d_one_field, get_sim_file_names('ra0h50dec-50', 90), '/home/nlharr/spt_code/bispec_routines/sptsz_input_scripts/test_generation/ra0h50dec-50_test_run_initial_90_apod_mask.sav','ra0h50dec-50', 90, 'test_trans2_store.fits'


pro trans_2d_one_field, sim_file, mask_file, field, band, output_file
sim_reso_arcmin = 0.25
map_reso_arcmin = 0.5
reso_rad = map_reso_arcmin*!dtor/60.

nsims = 100

bands = [90, 150, 220]
band_ind = where(bands eq band)

field_struct = get_field_struct()
fs = field_struct[where( field_struct.name eq field)]

;nx = fs.n1
;ny = fs.n2
;dnx = double(nx)
;dny = double(ny)

;loads the apodization mask
restore, mask_file
apmask = mask
maskfac_c = mean(apmask^2)

ap_size = size(apmask, /dimensions)
out_n0 = ap_size[0]
out_n1 = ap_size[1]

dout_n0 = double(out_n0)
dout_n1 = double(out_n1)

reso_rad = map_reso_arcmin*!dtor/60.

;loads the gridded cls and pls
theory_cls = get_highel2500_theory_dls(/cl)
ells = reform(theory_cls[0, *])

ells = [0,0,ells]
cls = reform(theory_cls[band_ind + 1, *])
cls = [0,0,cls]
pls = calc_pl(ells,sim_reso_arcmin)

print, 'inverting pls'
ipls = 1./pls

ellgrid=2.*!pi*make_fft_grid(reso_rad,out_n0,out_n1,fx=ugrid,fy=vgrid)
ellxgrid=2.*!pi*ugrid
ellygrid=2.*!pi*vgrid

thiscl2d = cls[ellgrid]
thisipl2d = ipls[ellgrid]

thispl2d = pls[ellgrid]

;loads the output simulations
simmaps_unpad = load_simmaps_from_file( sim_file, field)

simmaps_ps = dblarr(out_n0,out_n1,nsims)
sim_map = 0

;pad the sim map to be the same size as the apodization mask
simmaps = make_array(out_n0, out_n1, nsims, value = 0.0)
print, 'padding sim maps ', out_n0, ' ', out_n1
for k=0,nsims-1 do simmaps[*,*,k] = pad_2d_array(simmaps_unpad[*,*,k], out_n0, out_n1)

print, 'ffting sim maps'
for k=0,nsims-1 do simmaps_ps[*,*,k] = abs(fft(simmaps[*,*,k]*apmask))^2*reso_rad^2*dout_n0*dout_n1*1d12/maskfac_c
;        simmap_ps_tot = total(simmaps_ps,3)/100d0
tf_arr = fltarr(out_n0,out_n1,nsims)
whg = where(thiscl2d gt 0.)

print, 'formatting output'
for ii=0,nsims-1 do begin
    tftemp = fltarr(out_n0,out_n1)
    simmap_ps_tot = simmaps_ps[*,*,ii]
    tftemp[whg] = sqrt(simmap_ps_tot[whg]/thiscl2d[whg])
    tf_arr[*,*,ii] = tftemp * thispl2d ;scale by pixel window so we take account of it
endfor


tf_average = total(tf_arr,3);sum along sim axis
tf_average = tf_average/nsims


print, 'max_tf ', max(tf_average)
print, 'maskfac_c ', maskfac_c
print, 'dout_n0 ', dout_n0 
print, 'dout_n1 ', dout_n1 
print, 'reso_rad ', reso_rad 

;print, "WARNING NOT SCALING BY PIXEL WINDOW and NOT SAVING"

writefits,output_file,tf_average

end
