pro load_bispec_input_hdf5, in_file, map, reso_arcmin, apod_mask, winf_grid, ps_grid, cmb_ells, cmb_cls, field_name, bands, pt_clust_mask
  fid = H5F_OPEN(in_file)
  map = simple_load_hdf5( fid,  'map')
  reso_arcmin = simple_load_hdf5( fid,  'reso_arcmin')
reso_arcmin = reso_arcmin[0]
  apod_mask = simple_load_hdf5( fid,  'apod_mask')
  winf_grid = simple_load_hdf5( fid,  'winf_grid')
  ps_grid = simple_load_hdf5( fid,  'ps_grid')
  cmb_ells = simple_load_hdf5( fid,  'cmb_ells')
  cmb_cls = simple_load_hdf5( fid,  'cmb_cls')
  field_name = simple_load_hdf5( fid,  'field_name')
field_name = field_name[0]
  bands = simple_load_hdf5( fid,  'bands')
  pt_clust_mask = simple_load_hdf5( fid,  'pt_clust_mask')
  H5F_CLOSE,fid
end

