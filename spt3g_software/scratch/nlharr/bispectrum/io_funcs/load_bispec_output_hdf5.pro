pro load_bispec_output_hdf5, in_file, bispec, bispec_weight, field_name, bands, grid_bands, bispec_tag, ell_vec, delta_ell
  fid = H5F_OPEN(in_file)
  bispec = simple_load_hdf5( fid,  'bispec')
  bispec_weight = simple_load_hdf5( fid,  'bispec_weight')
  field_name = simple_load_hdf5( fid,  'field_name')
field_name = field_name[0]
  bands = simple_load_hdf5( fid,  'bands')
  grid_bands = simple_load_hdf5( fid,  'grid_bands')
  bispec_tag = simple_load_hdf5( fid,  'bispec_tag')
  ell_vec = simple_load_hdf5( fid,  'ell_vec')
  delta_ell = simple_load_hdf5( fid,  'delta_ell')
delta_ell = delta_ell[0]
  H5F_CLOSE,fid
end

