#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <fftw3.h>

#include <boost/shared_ptr.hpp>

#include "mapspectra/bispec_est.h"

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

/**
   map_pps:  map_pointer to pointers
   trip:  A bispectrum coordinate on the cube, contains the band and the el of the element
**/

//if you change this function you need to make sure to change the get_valid_trips function
//because that assumes that band is more important for setting the order.
int cmp_band_el(const void * x, const void * y ){
	band_el a = *((band_el*) x);
	band_el b = *((band_el*) y);
	if (a.band == b.band){
		if ( a.el < b.el ) return -1;
		else if ( a.el > b.el ) return 1;
		else  return 0;
	}
	else if (a.band > b.band) return 1;
	else return -1;
}


int cmp_trip(const void * x, const void * y ){
	trip a = *((trip*) x);
	trip b = *((trip*) y);
	int cmp = cmp_band_el( &(a.i[0]), &(b.i[0]));
	if (cmp == 0){
		cmp = cmp_band_el( &(a.i[1]), &(b.i[1]));
		if (cmp == 0) return cmp_band_el( &(a.i[2]), &(b.i[2]));
		else return cmp;
	}else return cmp;
}

trip create_trip_i(band_el a[3]){
	trip out;
	out.i[0] = a[0];
	out.i[1] = a[1];
	out.i[2] = a[2];
	qsort (&(out.i[0]), 3, sizeof(band_el), cmp_band_el);
	return out;
}



int sloppy_eq(double a, double b){
	const double EPS = 0.01;
	return (fabs(a-b) < EPS);
}


int is_tri_eq(trip t, double slop){
	double a[3];
	a[0] = t.i[0].el;
	a[1] = t.i[1].el;
	a[2] = t.i[2].el;
	return ((fabs(a[0]) + fabs(a[1]) + slop >= fabs(a[2])) &&
		(fabs(a[1]) + fabs(a[2]) + slop >= fabs(a[0])) &&
		(fabs(a[2]) + fabs(a[0]) + slop >= fabs(a[1])) );
}

int is_valid_trip(trip a, double slop){
	int is_sorted =  ((cmp_band_el(&(a.i[0]), &(a.i[1])) <= 0 ) && (cmp_band_el(&(a.i[1]), &(a.i[2])) <= 0   ));
	int tri_eq = is_tri_eq(a, slop);
	return is_sorted && tri_eq;
}

void print_int_arr(int * a, int num ){
	int i;
	for (i=0; i<num; i++) printf("a[%d]:%d\n",i,a[i]);
	printf("\n\n");
}

void print_trip_arr(trip * out_ts, int num ){
	int i;
	for (i=0; i<num; i++) printf( "Found trip: %f|%d %f|%d %f|%d\n", 
				      out_ts[i].i[0].el, out_ts[i].i[0].band, 
				      out_ts[i].i[1].el, out_ts[i].i[1].band, 
				      out_ts[i].i[2].el, out_ts[i].i[2].band );
	printf("\n\n");
}



//YOU NEED TO FREE THE RETURNED ARRAY
trip * get_valid_trips(int n_bands, int * bands, 
		       int n_els, double * els,
		       int * n_trips){
	int max_num = (n_bands*n_els)*(n_bands*n_els)*(n_bands*n_els);
	int bi,bj,bk,eli,elj,elk,i;
	
	
	double delta_el =  els[1]-els[0];
	for (i=2; i<n_els; i++){
		if (! sloppy_eq(els[i]-els[i-1], delta_el)){
			printf("delta el: %f, found delta %f\n",delta_el,els[i]-els[i-1]);
			printf("output el spacing needs to be the same, sorry.\n");
			exit(1);
		}
	}
	
	trip * trips;
	trips = (trip*) malloc( max_num * sizeof(trip) );
	
	int cur_trip = 0;
	
#ifdef DEBUG
	printf("generating full\n");
#endif
	
	//PREPARE FOR THE FOR LOOPS
	for (eli=0; eli< n_els; eli++)
		for (elj=eli; elj< n_els; elj++)
			for (elk=elj; elk< n_els; elk++)	
				
				for (bi=0; bi< n_bands; bi++)
					for (bj=0; bj< n_bands; bj++)
						for (bk=0; bk< n_bands; bk++){
							
							band_el be[3];
							be[0].band = bands[bi];
							be[0].band_ind = bi;
							be[0].el   = els[eli];
							be[0].el_ind = eli;
							
							be[1].band =bands[bj];
							be[1].band_ind = bj;
							be[1].el   = els[elj];
							be[1].el_ind = elj;
							
							be[2].band =bands[bk];
							be[2].band_ind = bk;
							be[2].el   = els[elk];
							be[2].el_ind = elk;
							
							trips[cur_trip] = create_trip_i( be );
							cur_trip ++;
						}
	
#ifdef DEBUG
	printf("filtering bad\n");
#endif
	
	// sets all the invalid trips to have the same triplet
	band_el bel;
	bel.band = -1; //being -1 is important because it means it will sort to the back
	bel.el   = -1.0f;
	for (i=0; i < cur_trip; i++){
		if ( !is_valid_trip( trips[i], delta_el)  ){
			trips[i].i[0] = bel;
		}
	}
	
	qsort (trips, cur_trip, sizeof(trip), cmp_trip);
	
	
	//flag redundant bolos
	for (i=0; i < cur_trip-1; i++){
		if ( cmp_trip(  &(trips[i]),&(trips[i+1]) ) == 0 ) trips[i].i[0] = bel;
	}
	
	//allocate memory for the returned list
	*n_trips = 0;
	for (i=0; i < cur_trip; i++){
		if (cmp_band_el(&(trips[i].i[0]),  &bel) != 0){
			*n_trips += 1;
		}
	}
	trip * out_trips;
	out_trips = (trip*) malloc(sizeof(trip)* (*n_trips));
	
	//build the returned list
	cur_trip = 0;
	for (i=0; i < *n_trips; i++){
		if (cmp_band_el(&(trips[i].i[0]),  &bel) != 0){
			out_trips[cur_trip] = trips[i];
			cur_trip++;
		}
	}
	
	
	//free the memory used while calculating
	free(trips);
	return out_trips;
	
}




bispec_t calc_bispec_est_sing_omp(map_t * m0, map_t * m1, map_t * m2, int n){
	if (n%4 != 0){
		printf("Due to loop unrolling your total map size needs to be a multiple of 4.  Just make your dimensions even.\n");
		exit(1);
	}
	
	bispec_t sum = 0.0;
	int i;
	for (i = 0; i < n; i += 4){
		sum += m0[i]*m1[i]*m2[i];
		sum += m0[i+1]*m1[i+1]*m2[i+1];
		sum += m0[i+2]*m1[i+2]*m2[i+2];
		sum += m0[i+3]*m1[i+3]*m2[i+3];
	}
	return sum;
}



//array of indices, [index*3 + first, sec, third]
//allocated bispectrum array
//pointer to pointers!
//map len


//requires out_bs be allocated and have enough memory for n_inds*sizeof(bispec_t)
//map_inds has length 3*n_inds since they are composed of triplets
void calc_bispec_est(  bispec_t * out_bs,
		       int * map_inds, int n_inds, 
		       map_t ** maps, int n_maps, 
		       int map_size){
	int i;
	//first check that the map_indices are valid
	for (i=0; i < n_inds * 3; i++){
		if (map_inds[i] >= n_maps){
			printf("invalid indices in calc_bispec_est code\n");
			printf("i:%d map_inds[i]:%d n_maps:%d\n", i, map_inds[i], n_maps);
			exit(1);
		}
	}
	
	//now do our thing
#ifdef OPENMP_FOUND
#pragma omp parallel for private(i) shared(out_bs, maps, map_size)
#endif
	for (i=0; i< n_inds; i++){
		out_bs[i] = calc_bispec_est_sing_omp(maps[map_inds[3*i]], 
						     maps[map_inds[3*i+1]], 
						     maps[map_inds[3*i+2]], map_size);
	} 
}


void print_sum_map(map_t * mp, int n){
	double sum = 0;
	int i;
	for(i=0; i< n; i++) sum +=mp[i];
	printf("found map sum %e\n", sum);
	
	
}

void print_sum_bispec(bispec_t * mp, int n){
	double sum = 0;
	int i;
	for(i=0; i< n; i++) sum +=mp[i];
	printf("found bsipec sum %e\n", sum);
	
}


//sometimes I like assinine names
//assumes bands are ordered
//you need to free map_pps
void fill_map_pps( map_t *** map_pps, int * n_bands,
		   map_t * ms0, map_t * ms1, map_t * ms2, 
		   int n_maps, int map_size, int bands[3] ){
	int i, j;

	//checks that the bands are sorted
	if (bands[0] > bands[1] || bands[1] > bands[2] ){
		printf("this code assumes the bands are sorted.\n");
		exit(1);
	}
  
	//find number of unique bands  
	*n_bands = 3;
	if (bands[0] == bands[1]) *n_bands -= 1;
	if (bands[1] == bands[2]) *n_bands -= 1;
  
	if (*n_bands == 2 && (bands[0] == bands[1])){
		printf("if you are passing two bands the code assumes they are 0 and 1\n");
		exit(1);
	}

	*map_pps = (map_t**) malloc(sizeof(map_t*) * (*n_bands) * n_maps);

	//printf("\n\nfound %d map_pps\n\n",  (*n_bands) * n_maps);

	//fills the map pps  
	map_t * mps[3];
	mps[0] = ms0;
	mps[1] = ms1;
	mps[2] = ms2;

	for (j=0; j < *n_bands; j++){
		for (i=0; i < n_maps; i++){
			//printf("\n\non j %d, i %d i*map_size: %d i+j*n_maps: %d\n\n", j,i, i*map_size, i+j*n_maps);
			(*map_pps)[i+j*n_maps] = &((mps[j])[i * map_size]);
		}
	}
}

//map_inds prealloced to 3*n_trips*sizeof(int)
void trip_to_mpps_ind(trip * trips, int n_trips, int * map_inds, int n_maps){
	int i;
	for (i=0; i < n_trips; i++){
		map_inds[i*3]   = trips[i].i[0].band_ind * n_maps + trips[i].i[0].el_ind;
		map_inds[i*3+1] = trips[i].i[1].band_ind * n_maps + trips[i].i[1].el_ind;
		map_inds[i*3+2] = trips[i].i[2].band_ind * n_maps + trips[i].i[2].el_ind;
	}
#ifdef DEBUG
	printf("n_maps: %d\n", n_maps);
#endif

}

//assumes all the pointers point to allocated memory of n_trips units
void trip_to_ret_inds(trip * trips, int n_trips,
		      int * band_0, int * band_1, int * band_2, 
		      int * el_0, int * el_1, int * el_2 ){
	int i;
	for (i=0; i < n_trips; i++){
		band_0[i] = trips[i].i[0].band;
		band_1[i] = trips[i].i[1].band;
		band_2[i] = trips[i].i[2].band;

		el_0[i] = trips[i].i[0].el_ind;
		el_1[i] = trips[i].i[1].el_ind;
		el_2[i] = trips[i].i[2].el_ind;
	}
}

int calc_full_bispec_est( boost::shared_ptr< std::vector<double> > ms0, 
			  boost::shared_ptr< std::vector<double> > ms1, 
			  boost::shared_ptr< std::vector<double> > ms2, 
			  int n_maps, int map_size,
			  std::vector<int> & bands, std::vector<double> &  els, 
			  //returned data, assumes n_ret alloced by calling, requires later freeing of other information
			  std::vector<double> & bispec_vals, //bispectrum values
			  std::vector<int> & band_0, 
			  std::vector<int> & band_1, 
			  std::vector<int> & band_2,  //bands
			  std::vector<int> & el_0, 
			  std::vector<int> & el_1, 
			  std::vector<int> & el_2 ){  //ells
	map_t ** map_pps = NULL;
	int n_bands;
	trip * trips;
	int n_trips;
	int * est_inds;

	//figures out how many unique maps there are an allocates/assigns the pointers to be in the 
	//appropriate index for use by the actual calculation code
	fill_map_pps( &map_pps, &n_bands,
		      &((*ms0)[0]), &((*ms1)[0]), &((*ms2)[0]), 
		      n_maps, map_size, &(bands[0]));
#ifdef DEBUG

	for (i=0; i< n_maps; i++) print_sum_map(map_pps[i], map_size);
	printf("getting valid trips \n");
	printf("n_bands:%d\n",n_bands);
	for(i=0; i<3; i++) printf("bands[%d]:%d\n",i, bands[i]);
  
#endif

	//gets the valid triplets
	trips =  get_valid_trips(n_bands, &(bands[0]), 
				 n_maps, &(els[0]), 
				 &n_trips);

#ifdef DEBUG
	printf("allocing memory\n");
#endif

	//allocates memory for the returned types 

	band_0 = std::vector<int>(n_trips, 0);
	band_1 = std::vector<int>(n_trips, 0);
	band_2 = std::vector<int>(n_trips, 0);
	el_0 = std::vector<int>(n_trips, 0);
	el_1 = std::vector<int>(n_trips, 0);
	el_2 = std::vector<int>(n_trips, 0);
	bispec_vals = std::vector<bispec_t>(n_trips, 0.0);
  
	//converts them to a form usable by the map est
#ifdef DEBUG
	printf("shuffling pointers\n");
#endif
  
	est_inds = (int *) malloc( sizeof(int) * 3 * n_trips);
	trip_to_mpps_ind(trips, n_trips,  est_inds, n_maps);
  
	//fills the returned indices
#ifdef DEBUG
	printf("generating returned indices\n");
#endif
	trip_to_ret_inds(trips, n_trips, &(band_0[0]), &(band_1[0]), &(band_2[0]), &(el_0[0]), &(el_1[0]), &(el_2[0]) );
  
	//actually calculates our estimator
#ifdef DEBUG
	printf("calc'ing estimator\n");
#endif
	calc_bispec_est(  &(bispec_vals[0]), est_inds, n_trips, map_pps, n_maps,  map_size);
  
  
#ifdef DEBUG
	printf("found %d ntrips\n", n_trips);
	print_sum_bispec(&(bispec_vals[0]), n_trips);
#endif

	free(est_inds);
	free(map_pps);
	free(trips);
	return 0;
}

//fmps   = get_filtered_ffts(fmp, els, delta_el, reso_rad)
//fmps = fft_maps(fmps, is_forward = False)
int fftw_is_init = 0;
void initialize_fftw(){
	if (!fftw_is_init){
		//fftw_init_threads();
		fftw_is_init = 1;
		//fftw_plan_with_nthreads(4);
	}
}

void finalize_fftw(){
	//fftw_cleanup_threads();
}

/**
//assumes that f_mps is allocated
//n0 is the slowest changing 
// n1 is fastest

f_mps output filtering maps
mp input map
map_scaling

map_scaling:  Scaling to apply to the map in frequency space (only applied if do_forward is true)
do_forward:  True if the input map is in map space

els:  the center of the ell rings to apply the fft (len of n_maps)
ell_grid:  The map of what ells the values are at each position in elspace
**/


void bispec_filter_maps(int n_maps, int n0, int n1, int do_forward,
			const std::vector<double> &  mp, 
			const std::vector<double> &  map_scaling, 
			const std::vector<double> &  ell_grid,
			const std::vector<double> &  els, double delta_el,
			std::vector<double> &  f_mps){
	initialize_fftw();
	fftw_complex * in;
	fftw_complex * out;
	fftw_complex * in_inv;
	fftw_complex * out_inv;

	int tmp_size = f_mps.size();
	if (tmp_size != n_maps * n0 * n1) f_mps = std::vector<double>(n_maps * n0 * n1, 0.0);

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0 * n1);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0 * n1);
	in_inv = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0 * n1);
	out_inv = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n0 * n1);

	fftw_plan forward_p;
	fftw_plan inverse_p;
	inverse_p = fftw_plan_dft_2d(n0, n1,in_inv, out_inv,
				     FFTW_BACKWARD, FFTW_ESTIMATE);
	if (do_forward){
		printf("doing forward\n");
		forward_p = fftw_plan_dft_2d(n0, n1, in, out,
					     FFTW_FORWARD, FFTW_ESTIMATE);
		for (int j=0; j< n0*n1; j++){
			in[j][0] = mp[j];
			in[j][1] = 0;
		}
		fftw_execute(forward_p);  
    
		for (int j=0; j< n0*n1; j++){
			out[j][0] *= map_scaling[j];
			out[j][1] *= map_scaling[j];
		}
	} else{
		for(int j=0; j< n0*n1; j++){ 
			out[j][0] = mp[j];
			out[j][1] = 0;
		}
	}

	for (int i=0; i < n_maps; i++){
		for (int j=0; j< n0*n1; j++){ 
			if (ell_grid[j] > els[i] - delta_el/2.0 && ell_grid[j] < els[i] + delta_el/2.0){
				in_inv[j][0] = out[j][0];
				in_inv[j][1] = out[j][1];
			}
			else{
				in_inv[j][0] = 0;
				in_inv[j][1] = 0;
			}
		}
		fftw_execute(inverse_p);
		for (int j=0; j< n0*n1; j++){
			f_mps[i*n0*n1 + j] = out_inv[j][0];
		}
	}

	fftw_free(in);
	fftw_free(in_inv);
	fftw_free(out);
	fftw_free(out_inv);
}

