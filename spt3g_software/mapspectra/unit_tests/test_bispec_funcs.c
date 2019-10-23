#include <stdio.h>
#include <string.h>
#include "CUnit/Basic.h"

#include "../c_bispec_est.h"
#include "../inpaint_quick.h"


/* The suite initialization function.
 * Opens the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */
int init_suite1(void)
{
  return 0;
}

/* The suite cleanup function.
 * Closes the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */
int clean_suite1(void)
{
  return 0;
}

/* Simple test of fprintf().
 * Writes test data to the temporary file and checks
 * whether the expected number of bytes were written.
 */
void testBandEl(void)
{

  band_el be[3];
  be[0].band   = 0;
  be[0].el     = 31.2f;
  be[0].el_ind = 0;

  be[1].band   = 0;
  be[1].el     = 32.2f;
  be[1].el_ind = 1;

  be[2].band   = 1;
  be[2].el     = 31.2f;
  be[2].el_ind = 0;

  CU_ASSERT( cmp_band_el( &(be[0]), &(be[0])) == 0 );
  CU_ASSERT( cmp_band_el( &(be[1]), &(be[0])) > 0 );  
  CU_ASSERT( cmp_band_el( &(be[0]), &(be[1])) < 0 );

  CU_ASSERT( cmp_band_el( &(be[2]), &(be[0])) > 0 );  
  CU_ASSERT( cmp_band_el( &(be[0]), &(be[2])) < 0 );

  CU_ASSERT( cmp_band_el( &(be[2]), &(be[1])) > 0 );  
  CU_ASSERT( cmp_band_el( &(be[1]), &(be[2])) < 0 );
}
void testTrips(void)
{

  band_el be[3];
  be[0].band   = 0;
  be[0].el     = 31.2f;
  be[0].el_ind = 0;

  be[2].band   = 0;
  be[2].el     = 32.2f;
  be[2].el_ind = 1;

  be[1].band   = 1;
  be[1].el     = 31.2f;
  be[1].el_ind = 0;
  

  band_el be2[3];
  be2[0].band   = 1;
  be2[0].el     = 31.2f;
  be2[0].el_ind = 0;

  be2[2].band   = 0;
  be2[2].el     = 32.2f;
  be2[2].el_ind = 1;

  be2[1].band   = 1;
  be2[1].el     = 31.2f;
  be2[1].el_ind = 0;
  


  band_el be3[3];
  be3[0].band   = 0;
  be3[0].el     = 31.2f;
  be3[0].el_ind = 0;

  be3[2].band   = 1;
  be3[2].el     = 32.2f;
  be3[2].el_ind = 1;

  be3[1].band   = 1;
  be3[1].el     = 31.2f;
  be3[1].el_ind = 0;

  band_el be4[3];
  be4[0].band   = 0;
  be4[0].el     = 31.2f;
  be4[0].el_ind = 0;

  be4[2].band   = 0;
  be4[2].el     = 32.2f;
  be4[2].el_ind = 1;

  be4[1].band   = 2;
  be4[1].el     = 31.2f;
  be4[1].el_ind = 0;


  band_el be5[3];
  be5[0].band   = 0;
  be5[0].el     = 2.5f;
  be5[0].el_ind = 0;

  be5[2].band   = 0;
  be5[2].el     = 2.51f;
  be5[2].el_ind = 0;

  be5[1].band   = 0;
  be5[1].el     = 5.0f;
  be5[1].el_ind = 1;



  trip t0 = create_trip_i(be);
  trip t1 = create_trip_i(be2);
  trip t2 = create_trip_i(be3);
  trip t3 = create_trip_i(be4);

  trip t4 = create_trip_i(be5);



  CU_ASSERT( is_tri_eq(t0, 0));
  CU_ASSERT( is_tri_eq(t1, 0));
  CU_ASSERT( is_tri_eq(t2, 0));
  CU_ASSERT( is_tri_eq(t3, 0));

  CU_ASSERT( is_tri_eq(t4, 0));

  CU_ASSERT( is_valid_trip(t0,0));
  CU_ASSERT( is_valid_trip(t1,0));
  CU_ASSERT( is_valid_trip(t2,0));
  CU_ASSERT( is_valid_trip(t3,0));

  CU_ASSERT( cmp_band_el(&(t0.i[0]), &(be[0])) == 0);
  CU_ASSERT( cmp_band_el(&(t0.i[1]), &(be[2])) == 0);
  CU_ASSERT( cmp_band_el(&(t0.i[2]), &(be[1])) == 0);

  CU_ASSERT( cmp_trip(  &t0, &t0) == 0);
  CU_ASSERT( cmp_trip(  &t1, &t0) > 0);
  CU_ASSERT( cmp_trip(  &t0, &t1) < 0);

  CU_ASSERT( cmp_trip(  &t2, &t0) > 0);
  CU_ASSERT( cmp_trip(  &t0, &t2) < 0);

  CU_ASSERT( cmp_trip(  &t3, &t0) > 0);
  CU_ASSERT( cmp_trip(  &t0, &t3) < 0);

  CU_ASSERT( cmp_trip(  &t1, &t2) > 0);
  CU_ASSERT( cmp_trip(  &t2, &t1) < 0);

  CU_ASSERT( cmp_trip(  &t2, &t3) > 0);
  CU_ASSERT( cmp_trip(  &t3, &t2) < 0);

}


void testGenTrips(void){


  size_t bands[2] = {0,1};

  size_t n_bands = 2;
  size_t n_els = 5;

  float els[] = { 2.5,  2.6, 2.7, 2.8,2.9,3.0,3.1, 3.2 };
  size_t n_trips;

  float els2[] = { 2.5, 7.5, 12.5 };

  float els3[] = { 2.5, 7.5, 12.5 };


  trip * out_ts;


  //expect to follow form n*(n-1)*(n-2)/6 + n*(n-1) + n
  //nbands then n_els
  out_ts =  get_valid_trips(1, bands, 2, els, &n_trips);
  CU_ASSERT( n_trips == 4);
  free(out_ts);


  out_ts =  get_valid_trips(1, bands, 8, els, &n_trips);
  CU_ASSERT( n_trips == 120);
  free(out_ts);

  out_ts =  get_valid_trips(2, bands, 8, els, &n_trips);
  CU_ASSERT( n_trips == 816);
  free(out_ts);

  /*
  out_ts =  get_valid_trips(1, bands, 3, els2, &n_trips);
  CU_ASSERT( n_trips == 6);
  free(out_ts);


  out_ts =  get_valid_trips(1, bands, 3, els3, &n_trips);
  CU_ASSERT( n_trips == 6);
  free(out_ts);
  */

  
  //qsort (out_ts, n_trips, sizeof(trip), cmp_trip);

  //printf("\n\nfound %zu trips\n\n", n_trips);
  //print_trip_arr(out_ts, n_trips );

  
}


void testCalcBiSing(void){
  map_t m0[] = {1, 2, 3, 5};
  map_t m1[] = {7, 11, 13, 17};
  map_t m2[] = {19, 23, 29, 31};

  size_t n = 4;

  CU_ASSERT(calc_bispec_est_sing_omp(m0, m1, m2, n) == 4405.);
  
}


void testCalcBiMult(void){


  map_t ms[12] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
  size_t n_maps = 3;
  size_t map_size = 4;

  map_t ** mpps;

  //mpps = (map_t**) malloc(n_maps * sizeof(map_t*));
  size_t n_bands;
  size_t bands[3] = {0,0,0};
  //fill_map_pps( mpps, ms, n_maps, map_size);

  fill_map_pps( &mpps, &n_bands,
		ms, NULL, NULL, 
		n_maps, map_size, bands );


  bispec_t out_bs[5];

  size_t map_inds[] = {0,1,2,  //4405
                       0,0,0,  //161
                       0,1,0, //593
                       1,0,0,
                       0,0,1};
  size_t n_inds = 5;

  printf("calcing bispec\n");
  calc_bispec_est(  out_bs,map_inds, n_inds, 
		    mpps, n_maps, 
		    map_size);  

  CU_ASSERT(out_bs[0] == 4405.);
  CU_ASSERT(out_bs[1] == 161.);
  CU_ASSERT(out_bs[2] == 593.);
  CU_ASSERT(out_bs[3] == 593.);
  CU_ASSERT(out_bs[4] == 593.);
 
 free( mpps);
  
}


void test_inpaint_algorithm(void){
  short mask_0[9] = { 0,1,1,
		      1,1,1,
		      1,1,0};

  float input_map_0[9] = { 12,0,0,
			   0,0,0,
			   0,0,48};

  float iter_1[9] = {12, 4, 0,
		     4, 0, 16,
		     0,16,48};



  float input_map_1[9] = { 12,0,0,
			   0,0,0,
			   0,0,48};



  float iter_2[9] = {12, 4, 10,
		     4, 10, 16,
		     10, 16, 48};

  int i,j;
  inpaint_map_laplace(mask_0, input_map_0, 3,3,1);
  for (i=0;i<9;i++){ 
    //printf("i: %d, expec: %f, actual; %f\n", i, iter_1[i], input_map_0[i]);
    CU_ASSERT_DOUBLE_EQUAL(iter_1[i],input_map_0[i], 0.01);
  }

  inpaint_map_laplace(mask_0, input_map_1, 3,3,2);
  for (i=0;i<9;i++){ 
    //printf("i: %d, expec: %f, actual; %f\n", i, iter_1[i], input_map_0[i]);
    CU_ASSERT_DOUBLE_EQUAL(iter_2[i],input_map_1[i], 0.01);
  }
 


} 




/* The main() function for setting up and running the tests.
 * Returns a CUE_SUCCESS on successful running, another
 * CUnit error code on failure.
 */
int main()
{
  CU_pSuite pSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS != CU_initialize_registry())
    return CU_get_error();

  /* add a suite to the registry */
  pSuite = CU_add_suite("Suite_Bispec", init_suite1, clean_suite1);
  if (NULL == pSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suite */
  /* NOTE - ORDER IS IMPORTANT - MUST TEST fread() AFTER fprintf() */
  if ((NULL == CU_add_test(pSuite, "test of band el ", testBandEl)) ||
      (NULL == CU_add_test(pSuite, "test of triplet funcs" , testTrips)) ||
      (NULL == CU_add_test(pSuite, "test of triplet generation" , testGenTrips)) ||
      (NULL == CU_add_test(pSuite, "test of estimator on multiple bits" , testCalcBiMult)) ||
      (NULL == CU_add_test(pSuite, "test of estimator ", testCalcBiSing)) ||
      (NULL == CU_add_test(pSuite, "test of inpainting alg ", test_inpaint_algorithm))
      )
    {
      CU_cleanup_registry();
      return CU_get_error();
    }

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}
