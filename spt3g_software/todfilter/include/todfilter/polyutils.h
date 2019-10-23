#pragma once
#include <vector>
#include <G3Frame.h>
#include <G3Module.h>
#include <G3Timestream.h>
#include <G3Map.h>
#include <serialization.h>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <coordinateutils/G3SkyMap.h>

typedef std::shared_ptr<gsl_matrix> GslMatrixPtr;
typedef std::shared_ptr<const gsl_matrix> GslMatrixConstPtr;
GslMatrixPtr get_shared_gsl_matrix(size_t n1, size_t n2);

typedef std::shared_ptr<gsl_vector> GslVectorPtr;
typedef std::shared_ptr<const gsl_vector> GslVectorConstPtr;
GslVectorPtr get_shared_gsl_vector(size_t n);


void get_thin_svd_of_matrix(const GslMatrixPtr A,
			    GslMatrixPtr & U, 
			    GslMatrixPtr & V, 
			    GslVectorPtr & S);

void filter_from_model_svd( GslVectorPtr y_in, const GslMatrixPtr mod,
			    const GslMatrixPtr U, 
			    const GslMatrixPtr V, const GslVectorPtr S);



struct FilterMask : public G3FrameObject{
	G3MapInt has_masked_pixels;
	G3MapVectorInt pixel_mask;
	std::string Description() const{return std::string("FilterMask");}
	template <class A> void serialize(A &ar, unsigned v){
		G3_CHECK_VERSION(v);
		ar & cereal::make_nvp("base", cereal::base_class<G3FrameObject>(this));
		ar & cereal::make_nvp("has_masked_pixels", has_masked_pixels);
		ar & cereal::make_nvp("pixel_mask", pixel_mask);
	}
	bool is_valid() const{
		for (auto it = pixel_mask.begin(); it != pixel_mask.end(); it++){
			bool pmask_val = has_masked_pixels.at(it->first);
			if (pmask_val){
				int mask_val = 0;
				for (size_t i = 0; i < it->second.size(); i++) {
					mask_val |= (it->second.at(i));
				}
				if (!mask_val) {
					return false;
				}
			}
		}
		return true;
	}

};
G3_POINTERS(FilterMask);
G3_SERIALIZABLE(FilterMask, 1);


///////////////////////////////////
///////////////////////////////////


FilterMaskPtr make_empty_filter_mask(G3TimestreamMapConstPtr gm);
FilterMaskPtr make_filter_mask(G3SkyMapConstPtr point_source_mask,
			       G3MapVectorIntConstPtr detector_pointing);


G3TimestreamMapPtr poly_and_mhpf_filter_ts_data(G3TimestreamMapConstPtr timestreams, 
						FilterMaskConstPtr mask, 
						double sample_frequency, double freq_cutoff, 
						int poly_order);


G3TimestreamMapPtr notch_filter_lls(
    G3TimestreamMapConstPtr timestreams,
    FilterMaskConstPtr mask,
    const std::vector<double> & freqs);


G3TimestreamMapPtr poly_filter_ts_data_with_abscissa(
	G3TimestreamConstPtr abscissa,
	G3TimestreamMapConstPtr timestreams, 
	FilterMaskConstPtr mask, 
	int poly_order,	int index_poly_order);


void rolling_mean_filter_ts_data(G3TimestreamMapConstPtr timestreams, int filter_width,
				 G3TimestreamMapPtr out_timestreams);


