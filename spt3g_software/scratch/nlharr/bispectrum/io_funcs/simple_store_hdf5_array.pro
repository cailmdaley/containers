
pro simple_store_hdf5_array, fid, data, data_name
  print, 'simple_store_hdf5_array ',data_name
  print, size(data, /DIMENSION)
  datatype_id = H5T_IDL_CREATE(data)
  dataspace_id = H5S_CREATE_SIMPLE(size(data,/DIMENSIONS))
  ;; create dataset in the output file
  dataset_id = H5D_CREATE(fid,data_name,datatype_id,dataspace_id)
  ;; write data to dataset
  H5D_WRITE,dataset_id,data
  ;; close all open identifiers
  H5D_CLOSE,dataset_id
  H5S_CLOSE,dataspace_id
  H5T_CLOSE,datatype_id
end
