

;test command;
;make_complete_one_field, 'ra5h30dec-55','/home/nlharr/spt_code/bispec_routines/sptsz_input_scripts/test_generation/','test_run_initial', 90, 5, 5, 3, 1, 2.5e-14

pro make_complete_one_field, field, out_dir, tmp_tag_ini, band, $
                             cluster_mask_size, big_ptsrc_mask_size, small_ptsrc_mask_size, $
                             map_size_increase, cluster_mass_cut, no_update


;the mask size is in units of pixels away from the center we mask, 
;where the size of the pixels we mask is 1 + 2*size.  So size = 2
;means masking a five pixel square 
;  mask sizes is given in pixels
;
;
;
; apod mask type
; inpainting type  ; how to vary it
; mask sizes
; 






print, 'map_size_inc', map_size_increase


reso_arcmin = 0.5

field_struct = get_field_struct()
fs = field_struct[where( field_struct.name eq field)]

field_obs_str = get_field_obs_lst()

;for the bands 0=90, 1 = 150, 2 = 220


if band eq 90 then begin 
    band_ind = 0
endif else if band eq 150 then begin 
    band_ind = 1
endif else if band eq 220 then begin
    band_ind = 2
endif else begin
    print, "pick a reasonable band"
    return
end


tmp_tag = tmp_tag_ini + '_'+strtrim(band,2)


sim_file_name = get_sim_file_names(field, band)


coadd_fn_out = out_dir+'/'+field+'_'+tmp_tag+'_coadd.sav'
apod_fn_out = out_dir+'/'+field+'_'+tmp_tag+'_apod_mask.sav'
dummy_noise_fn_out = out_dir+'/'+field+'_'+tmp_tag+'_dummy_noise.fits'
noise_fn_out = out_dir+'/'+field+'_'+tmp_tag+'_noise_psd.fits'

trans2d_fn_out = out_dir+'/'+field+'_'+tmp_tag+'_trans2d.fits'


out_hdf5_fn_out = out_dir+'/'+field+'_'+tmp_tag+'_bundle.hdf5'
;map, apod, noise psd,

;coadd_one_field, field, 
ptsrc_file_name = '/home/nlharr/Bispec_2500/pt_src_lists/ptsrc_config_'+field+'.txt'
cluster_file = '/home/nlharr/Bispec_2500/pt_src_lists/spt_clust_allz_list_07mar2014.txt'
cluster_mass_cut = 2.5e-14





print, 'ptsrcing', ptsrc_file_name
make_full_mask_lst, ptsrc_file_name, field, cluster_file, cluster_mass_cut, full_ra, full_dec, full_desc
pt_clust_mask = convert_mask_lst_to_mask( full_ra, full_dec, full_desc, field, cluster_mask_size, big_ptsrc_mask_size, small_ptsrc_mask_size)




if ~ keyword_set(no_update) then no_update = 0

if (~ file_test(coadd_fn_out)) && ~ no_update then begin
    print, 'coadding'
    coadd = spt_coadd_dir('/data23/lizinvt/highel2500/coadds/combmapsx2_jittertimeordered_20130616_good/'+field+'/', 'map_'+field+'_'+strtrim(band,2)+'*')
    save, coadd, filename = coadd_fn_out
endif

;get the apodization mask
if (~ file_test(apod_fn_out)) && ~ no_update then begin
    print, 'apoding'
    apod_mask = apod_one_field( field, coadd_fn_out, apod_fn_out, map_size_increase)
endif
;generate the noise
print, 'noising'
if (~ file_test(noise_fn_out)) && ~ no_update then begin
    noise_one_field, field, band, apod_fn_out, noise_fn_out, dummy_noise_fn_out
endif
;trans 2d
print, 'trans_2d'
if (~ file_test(trans2d_fn_out)) && ~ no_update then begin
    trans_2d_one_field, sim_file_name, apod_fn_out, field, band, trans2d_fn_out
endif

;get the point source mask


print, 'getting cl spec'
;reads in the cls
cl_spec = get_highel2500_theory_dls(/cl)


;read in the apodization mask
restore, apod_fn_out
apod_mask = mask
apod_dims = size(apod_mask, /dimension)

;reading in the map
;fts = read_spt_fits(coadd_fn_out)
;fts = spt_expand_smallmap(fts)

restore, coadd_fn_out
out_map = pad_2d_array(coadd.map, apod_dims[0], apod_dims[1]) * 1e6 ; pads the map and scales to uK units

pt_clust_mask = pad_2d_array(pt_clust_mask, apod_dims[0], apod_dims[1]) ; pads the map and scales to uK units

print, 'pt_clust_mask max ', max(pt_clust_mask), ' min ', min(pt_clust_mask)


;loads the winf_grid
winf_grid = readfits(trans2d_fn_out)

;loads the ps_grid
restore, noise_fn_out
ps_grid = psd*1e6
;pro store_bispec_input_hdf5, out_file, map, reso_arcmin, apod_mask, winf_grid, ps_grid, cmb_ells, cmb_cls, field_name, bands, pt_clust_mask


store_bispec_input_hdf5, out_hdf5_fn_out, out_map, reso_arcmin, apod_mask, winf_grid, ps_grid, cl_spec[0,*], cl_spec[band_ind+1,*], field, [band], pt_clust_mask
end
