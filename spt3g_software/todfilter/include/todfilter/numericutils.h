#pragma once

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <G3Module.h>
#include <G3Map.h>
#include <G3Timestream.h>
/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 *  with modifications by nlharr to make it use templates
 */



G3Timestream mean_filter_timestream(const G3Timestream & ts);
double get_min(const G3Timestream & ts);
double get_max(const G3Timestream & ts);
void sum_reduce_vec_to_size(double * iv, size_t iv_size, size_t new_size);


template <typename T>
void elem_swap(T & a, T & b){
  T t=(a);(a)=(b);(b)=t;
}

//returns the lower median of the group
//modifies the input array
template <typename elem_type>
elem_type quick_select(elem_type arr[], int n) 
{
  int low, high ;
  int median;
  int middle, ll, hh;

  low = 0 ; high = n-1 ; median = (low + high) / 2;
  for (;;) {
    if (high <= low) /* One element only */
      return arr[median] ;

    if (high == low + 1) {  /* Two elements only */
      if (arr[low] > arr[high])
	elem_swap<elem_type>(arr[low], arr[high]) ;
      return arr[median];
    }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    elem_swap<elem_type>(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       elem_swap<elem_type>(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     elem_swap<elem_type>(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    elem_swap<elem_type>(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do ll++; while (arr[low] > arr[ll]) ;
      do hh--; while (arr[hh]  > arr[low]) ;

      if (hh < ll)
        break;

      elem_swap<elem_type>(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    elem_swap<elem_type>(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
      low = ll;
    if (hh >= median)
      high = hh - 1;
  }
}


//Computes the median absolute deviation if that wasn't fucking obvious
template <typename T>
T compute_median_absolute_deviation(const T arr[], int n){
  T * new_arr = (T*) malloc(sizeof(T) * n);
  if (! new_arr){
    printf("malloc error in compute_median_absolute_deviation\n");
    exit(1);
  }
  memcpy(new_arr, arr, n * sizeof(T));
  T median = quick_select(new_arr, n);
  memcpy(new_arr, arr, n * sizeof(T));
  for (int i=0; i < n; i++) new_arr[i] = fabs(new_arr[i] -  median);
  T mad = quick_select(new_arr, n);
  free(new_arr);
  return mad;  
}

double get_mad_to_variance_correction_normal_dist();



template <typename T>
void find_local_maxima(const std::vector<T> & in_vec, 
		       int min_index, size_t max_index,
		       int & out_index, T & out_val){
  min_index = min_index < 0 ? 0 : min_index;
  max_index = max_index >= in_vec.size() ? in_vec.size() - 1 : max_index;

  out_index = max_index;
  out_val = in_vec[max_index];
  for (size_t i=min_index; i < max_index; i++){
    if ( in_vec[i] > out_val){
      out_val = in_vec[i];
      out_index = i;
    }
  }
}

#define MAD_VARIANCE_ADDER_DOCSTR \
	"Uses the median absolute deviation to estimate the variance of each timestream " \
        "and adds it to the frame.  This assumes gaussian noise.\n\n" \
	"ts_key: [->G3TimestreamMap] timestream key\n\n" \
	"variance_output_key [->G3MapDouble] where the variance is stored\n\n" 

class MadVarianceAdder : public G3Module{
 public:
  MadVarianceAdder(std::string ts_key,std::string variance_output_key);
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string variance_output_key_;
  std::string ts_key_;
};

#define VARIANCE_ADDER_DOCSTR \
        "Estimates the variance of each timestream " \
        "and adds it to the frame. \n\n" \
        "ts_key: [->G3TimestreamMap] timestream key\n\n" \
        "variance_output_key [->G3MapDouble] where the variance is stored\n\n"

class VarianceAdder : public G3Module{
 public:
  VarianceAdder(std::string ts_key,std::string variance_output_key);
  void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
 private:
  std::string variance_output_key_;
  std::string ts_key_;
};


