pro store_bispec_input_hdf5, out_file, map, reso_arcmin, apod_mask, winf_grid, ps_grid, cmb_ells, cmb_cls, field_name, bands, pt_clust_mask
  fid = H5F_CREATE(out_file)
  simple_store_hdf5_array, fid, map, 'map'
  simple_store_hdf5_scalar, fid, reso_arcmin, 'reso_arcmin'
  simple_store_hdf5_array, fid, apod_mask, 'apod_mask'
  simple_store_hdf5_array, fid, winf_grid, 'winf_grid'
  simple_store_hdf5_array, fid, ps_grid, 'ps_grid'
  simple_store_hdf5_array, fid, cmb_ells, 'cmb_ells'
  simple_store_hdf5_array, fid, cmb_cls, 'cmb_cls'
  simple_store_hdf5_scalar, fid, field_name, 'field_name'
  simple_store_hdf5_array, fid, bands, 'bands'
  simple_store_hdf5_array, fid, pt_clust_mask, 'pt_clust_mask'
  H5F_CLOSE,fid
end

