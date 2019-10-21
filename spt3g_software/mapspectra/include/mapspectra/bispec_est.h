#ifndef C_BISPEC_EST_H
#define C_BISPEC_EST_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>

//the types calculated are defined here so it is easy to change in the future
typedef double map_t;
typedef double bispec_t;

struct band_el{
  int band;
  int band_ind;
  float el;
  int el_ind;
};

typedef struct band_el band_el;//because this comes from c

struct trip{
  band_el i[3];
};

typedef struct trip trip; //because this comes from c


/**
The only two public functions

int calc_full_bispec_est(//input data
			  map_t * ms0, map_t * ms1, map_t * ms2,
			  int n_maps, int map_size,
			  int * bands, float * els,
			  //returned data
			  int * n_ret,
			  bispec_t ** bispec_vals,
			  int ** band_0, int ** band_1, int ** band_2, 
			  int ** el_0, int ** el_1, int ** el_2 );
 **/

void bispec_filter_maps(int n_maps, int n0, int n1, int do_forward,
			const std::vector<double> &  mp, 
			const std::vector<double> &  map_scaling, 
			const std::vector<double> &  ell_grid,
			const std::vector<double> &  els, double delta_el,
			std::vector<double> &  f_mps);


////////////////////////////////////////
////////////////////////////////////////
///Internal functions///////////////////
/////Not for public consumption/////////
////////////////////////////////////////


//cmp functions return positive if x>y, negative if x<y, 0 if equal
int cmp_band_el(const void * x, const void * y );
int cmp_trip(const void * x, const void * y );
trip create_trip_i(band_el a[3]);
int is_tri_eq(trip t, float slop);
int is_valid_trip(trip a, float slop);
trip * get_valid_trips(int n_bands, int * bands, 
		       int n_els, float * els,
		       int * n_trips);
bispec_t calc_bispec_est_sing_omp(map_t * m0, map_t * m1, map_t * m2, int n);
void calc_bispec_est(  bispec_t * out_bs, int * map_inds, int n_inds, 
		       map_t ** maps, int n_maps, int map_size);
void fill_map_pps( map_t *** map_pps, int * n_bands,
		   map_t * ms0, map_t * ms1, map_t * ms2, 
		   int n_maps, int map_size, int bands[3] );
void print_trip_arr(trip * out_ts, int num );


#endif
