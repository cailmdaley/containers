pro store_bispec_field_output_hdf5, out_file, bispec, bispec_variance, bispec_weight_used, field_name, bands, grid_bands, bispec_tag, ell_vec, delta_ell
  fid = H5F_CREATE(out_file)
  simple_store_hdf5_array, fid, bispec, 'bispec'
  simple_store_hdf5_array, fid, bispec_variance, 'bispec_variance'
  simple_store_hdf5_array, fid, bispec_weight_used, 'bispec_weight_used'
  simple_store_hdf5_scalar, fid, field_name, 'field_name'
  simple_store_hdf5_array, fid, bands, 'bands'
  simple_store_hdf5_array, fid, grid_bands, 'grid_bands'
  simple_store_hdf5_array, fid, bispec_tag, 'bispec_tag'
  simple_store_hdf5_array, fid, ell_vec, 'ell_vec'
  simple_store_hdf5_scalar, fid, delta_ell, 'delta_ell'
  H5F_CLOSE,fid
end

