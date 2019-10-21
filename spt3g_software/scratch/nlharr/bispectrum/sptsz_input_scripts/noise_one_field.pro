


pro noise_one_field, field_name, band, apod_mask_file, out_file, jack_file
reso_arcmin = 0.5
njack = 100
fits_map_names = get_coadd_file_list(field_name, band)
make_psds_wrapper,'',maplist=fits_map_names, band, apod_mask_file, njack, out_file,reso_arcmin=reso_arcmin,jackfile=jack_file
end
