

pro coadd_one_field, field_name, band, file_out

reso_arcmin = 0.5
fits_map_names = get_coadd_file_list(field_name, band)

coadd_fits_maps, filelist=fits_map_names, fileout=file_out, band=band, /uniform

end
