
function simple_load_hdf5, fid,  d_name
  dataset_id = H5D_OPEN(fid, d_name)
  val = H5D_READ(dataset_id)
  H5D_CLOSE, dataset_id
  return, val
end
