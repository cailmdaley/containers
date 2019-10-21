#include <pybindings.h>
#include <stdlib.h>
#include <vector>
#include <string> 

#include <G3.h>
#include <G3Frame.h>

#include <hdf5.h>
#include <hdf5_hl.h>




hid_t open_hdf5_read(std::string fname){
  return H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
}

void close_hdf5_file(hid_t hdf5_file){
  H5Fclose( hdf5_file ) ;
}

std::vector<uint8_t> read_hdf5_bitfield(hid_t hdf5_file, std::string path,
			       std::vector<int> & dim_arr ){
  //if (H5LTpath_valid ( hdf5_file, path.c_str(), 1)){
  if(1){
    hid_t dset, space;
    herr_t status;
    int ndims, k;
    size_t nelems;
    hsize_t * temp_dim_arr;
    dset = H5Dopen (hdf5_file, path.c_str(), H5P_DEFAULT);
    space = H5Dget_space(dset); 
    //check the dimensions we think it is are the output dimensions so we don't 
    //have some overflow
    ndims = H5Sget_simple_extent_ndims(space);


    temp_dim_arr = (hsize_t *) malloc( sizeof(hsize_t) * ndims);
    dim_arr = std::vector<int>(ndims > 0 ? ndims : 1);
    
    //loads the dimensions in
    H5Sget_simple_extent_dims(space, temp_dim_arr, NULL );
    for (k=0; k < ndims; k++) dim_arr[k] = (size_t)temp_dim_arr[k];
    if (ndims==0)  dim_arr[0] = 1;
    
    
    nelems = H5Sget_simple_extent_npoints( space );
    std::vector<uint8_t> buffer( nelems * sizeof(uint8_t));
    //reads the data
    status = H5Dread( dset, H5T_NATIVE_B8, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(buffer[0]) );
    //free the memory
    status = H5Dclose(dset);
    status = H5Sclose(space);    
    free(temp_dim_arr);

    return buffer;
  }else{
    return std::vector<uint8_t>();
  }
}



namespace bp = boost::python;
PYBINDINGS("util"){
  bp::def("open_hdf5_read", open_hdf5_read);
  bp::def("close_hdf5_file", close_hdf5_file);

  bp::def("read_hdf5_bitfield", read_hdf5_bitfield, 
	  "This is a really hack solution to read bitfields from spt idfs.  "
	  "Will probably cause all the problems in the future.");
}
