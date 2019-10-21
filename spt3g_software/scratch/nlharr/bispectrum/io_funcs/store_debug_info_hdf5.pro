pro store_debug_info_hdf5, out_file, itfgrid, fftm_real, fftm_imag, ellgrid
  fid = H5F_CREATE(out_file)
  simple_store_hdf5_array, fid, itfgrid, 'itfgrid'
  simple_store_hdf5_array, fid, fftm_real, 'fftm_real'
  simple_store_hdf5_array, fid, fftm_imag, 'fftm_imag'
  simple_store_hdf5_array, fid, ellgrid, 'ellgrid'
  H5F_CLOSE,fid
end

