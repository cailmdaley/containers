#include <pybindings.h>
#include <serialization.h>
#include <todfilter/polyutils.h>
#include <todfilter/numericutils.h>
#include <G3Map.h>

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <iostream>
#include <vector>
#include <complex>


#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


#ifdef OPENMP_FOUND
#include <omp.h>

#define FFTW_OMP
#define LLS_OMP
#endif

#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))


G3_SERIALIZABLE_CODE(FilterMask);

//lapack is not currently set to compile. If you want to use lapack you'll need to modify
//the build system
//#define ACTUALLY_USE_LAPACK 

//sptsz uses       lpfilt[i]=exp(-1*pow(freqs[i]/lowpass,6));
//And since when did we ever try new things?


//////////////////////////////////////////////////
/////// Linear Least Squares Filter //////////////
//////////////////////////////////////////////////
//////And useless comments.../////////////////////
//////////////////////////////////////////////////

GslMatrixPtr get_shared_gsl_matrix(size_t n1, size_t n2){ 
	return GslMatrixPtr(gsl_matrix_calloc (n1, n2), gsl_matrix_free);
}

GslVectorPtr get_shared_gsl_vector(size_t n){ 
	return GslVectorPtr(gsl_vector_calloc(n), gsl_vector_free);
}

void matrix_info( gsl_matrix * m, const char *info){
	log_debug("Got Matrix %s : M:%zu N:%zu\n", info, m->size1, m->size2);
}

void get_thin_svd_of_matrix(const GslMatrixPtr A,
			    GslMatrixPtr & U, 
			    GslMatrixPtr & V, 
			    GslVectorPtr & S){
	log_debug("Enter"); 
	size_t M = A->size1;
	size_t N = A->size2;
	gsl_vector * work = gsl_vector_alloc(N);
	U = get_shared_gsl_matrix(M,N);
	gsl_matrix_memcpy(U.get(), A.get());
	V = get_shared_gsl_matrix(N,N);
	S = get_shared_gsl_vector(N);
	gsl_linalg_SV_decomp (U.get(), V.get(), S.get(), work) ;  //A is now U
	gsl_vector_free(work);
	log_debug("Exit");
}

void get_qr_of_matrix(const GslMatrixPtr A,
		      GslMatrixPtr & qr, 
		      GslVectorPtr & tau){
	log_debug("Enter"); 
	size_t M = A->size1;
	size_t N = A->size2;
	tau = get_shared_gsl_vector(N);
	qr = get_shared_gsl_matrix(M,N);
	gsl_matrix_memcpy(qr.get(), A.get());
	gsl_linalg_QR_decomp (qr.get(), tau.get());
	log_debug("Exit");
}


//stubs the filter mask
FilterMaskPtr make_empty_filter_mask(G3TimestreamMapConstPtr gm){
	FilterMaskPtr om = boost::make_shared<FilterMask>();
	for (auto iter = gm->begin(); iter != gm->end(); iter++){
		om->has_masked_pixels[iter->first] = 0;
	}
	return om;
}

//For filter masks, 0 is not masked, 1 is masked
//makes a filter mask from detector pointing
FilterMaskPtr make_filter_mask(G3SkyMapConstPtr point_source_mask,
			       G3MapVectorIntConstPtr detector_pointing){
	log_debug("making filter mask");
	FilterMaskPtr out_mask(new FilterMask);
	for (auto iter=detector_pointing->begin(); iter!=detector_pointing->end(); iter++){
		out_mask->pixel_mask[iter->first] = 
			std::vector<int>(detector_pointing->at(iter->first).size());
		std::vector<int> & dvec = out_mask->pixel_mask[iter->first];
		const std::vector<int> & pvec = detector_pointing->at(iter->first);
		int psm_sum = 0;
		for (unsigned int i = 0; i < pvec.size(); i++){
			int val =  (*point_source_mask)[pvec[i]];
			dvec[i] = val;
			psm_sum += val;
		}
		out_mask->has_masked_pixels[iter->first] = ! ! psm_sum;
	}
	return out_mask;
}

//makes the masked model matrix from a model matrix and an in matrix
void make_mask_lls_model_matrix( const GslMatrixPtr in_mat,
				 const std::vector<int> & mask,
				 unsigned int n_params,
				 GslMatrixPtr & out_mat){
	out_mat = get_shared_gsl_matrix(in_mat->size1, in_mat->size2  );
	memcpy(out_mat->data, in_mat->data, in_mat->size1 * in_mat->size2 * sizeof(double));
	for (unsigned int i=0; i < mask.size(); i++){
		if (mask[i]){
			for (unsigned int j=0; j < n_params; j++){
				out_mat->data[j + i * n_params] = 0;
			}
		}
	}
}

GslVectorPtr get_fit_params(double * data, 
			    size_t len,
			    GslMatrixConstPtr U, 
			    GslMatrixConstPtr V, 
			    GslVectorConstPtr S) {
	gsl_vector * gsl_v = gsl_vector_alloc(len);
	gsl_v->size = len;
	gsl_v->stride = 1;
	memcpy(gsl_v->data, data, sizeof(double) * len);
	gsl_v->block = NULL;
	gsl_v->owner = 0;
	
	size_t M = U->size1;
	size_t N = U->size2;
	
	GslVectorPtr x = get_shared_gsl_vector(N);
	int status = gsl_linalg_SV_solve (U.get(), V.get(), S.get(), gsl_v, x.get());
	gsl_vector_free(gsl_v);
	return x;
}



GslMatrixPtr make_polarization_fit( G3SkyMapPtr Q,
				    G3SkyMapPtr U,
				    G3SkyMapPtr mask
	) {
	size_t len = Q->size();
	GslMatrixPtr m = get_shared_gsl_matrix( len, 2);
	//memcpy(m->data + 0 * len, &((*Q)[0]), sizeof(double) * len);
	//memcpy(m->data + 1 * len, &((*U)[0]), sizeof(double) * len);
	for (size_t j=0; j < len; j++) {
		if (!((*mask)[j])) {
			for (size_t i=0; i < 2; i++) {
				m->data[i + j*2] = 0;
			}
		} else {
			m->data[0 + j*2] = (*Q)[j];
			m->data[1 + j*2] = (*U)[j];
		}
	}
	return m;
}

G3VectorDouble get_polarization_fit_params(G3SkyMapPtr Q,
					   G3SkyMapPtr U,
					   G3SkyMapPtr mask,
					   G3SkyMapPtr map) {

	GslMatrixPtr fit_mat = make_polarization_fit( Q,U,mask);
	GslMatrixPtr U_m; 
	GslMatrixPtr V_m; 
	GslVectorPtr S_m;
	get_thin_svd_of_matrix(fit_mat,U_m,V_m,S_m);
	GslVectorPtr fit_params = get_fit_params( &((*map)[0]), map->size(),
						  U_m,V_m,S_m);
	G3VectorDouble out(2,0);
	out[0] = fit_params->data[0];
	out[1] = fit_params->data[1];
	return out;
}


void simple_common_mode_filter_cpp(G3TimestreamMapPtr ts_map,
				   const std::vector<std::string> & ts_map_keys){
	std::list<std::string> keys;
	for (auto it=ts_map_keys.begin(); it != ts_map_keys.end(); it++) {
		if ( ts_map->find(*it) != ts_map->end() ) {
			keys.push_back(*it);
		}
	}
	if (keys.size() > 0){
		G3Timestream ts_sum(*(ts_map->begin()->second));
		ts_sum *= 0;
		for (auto it=keys.begin(); it!=keys.end(); it++){
			ts_sum += *((*ts_map)[*it]);
		}
		ts_sum /= (double) keys.size();
		for (auto it=keys.begin(); it!=keys.end(); it++){
			(*(*ts_map)[*it]) -= ts_sum;
		}
	}
}


void masked_common_mode_filter_cpp(G3TimestreamMapPtr ts_map,
				   const G3TimestreamMapPtr ts_map_masked,
				   const std::vector<std::string> & ts_map_keys){
	std::list<std::string> keys;
	for (auto it=ts_map_keys.begin(); it != ts_map_keys.end(); it++) {
		if ( ts_map->find(*it) != ts_map->end() ) {
			keys.push_back(*it);
		}
	}
	if (keys.size() > 0){
		G3Timestream ts_sum(*(ts_map->begin()->second));
		ts_sum *= 0;
		for (auto it=keys.begin(); it!=keys.end(); it++){
			ts_sum += *((*ts_map_masked)[*it]);
		}
		ts_sum /= keys.size();
		for (auto it=keys.begin(); it!=keys.end(); it++){
			(*(*ts_map)[*it]) -= ts_sum;
		}
	}
}


GslMatrixPtr make_poly_and_highpass_fit( double sample_frequency, double freq_cutoff, 
					 int poly_order, int len){
	g3_assert(poly_order >= 0 || freq_cutoff > 0);
	
	poly_order = poly_order < 0 ? 0 : poly_order;
	int nharms =  floor(freq_cutoff / sample_frequency * len);
	if (nharms < 0) nharms = 0;
	
	double effdf = 2.0 * M_PI  / ((double) len);
	
	int n_params = 2 * nharms + 1 + poly_order;
	GslMatrixPtr m = get_shared_gsl_matrix( len, n_params);
	double * out_model_matrix = m->data;
	
	//poly or hpf always have a mean filter
	for (int i=0; i < len; i++){
		out_model_matrix[n_params * i] = 1.0f;
	} 
	//cover the first case
	if (poly_order >0){
		int j = 1;
		for (int i=0; i < len; i++){
			double x = -1.0f + 2 * ((double)i)/((double) len - 1.0f);
			out_model_matrix[n_params * i + j] = x;
		}
	}
	
	for (int j = 2; j <  1 + poly_order; j++){
		for (int i=0; i < len; i++){
			double x = (-1.0f + 2 * ((double)i)/((double) len - 1.0f));
			double j_leg = j;
			out_model_matrix[n_params * i + j] = //pow(x,j);
				((2*j_leg-1) * x * out_model_matrix[(j-1)+i*n_params] -
				 (j_leg-1) * out_model_matrix[(j-2)+i*n_params])/j_leg;
		}
	}
	for (int j = poly_order + 1; j < 2*nharms + poly_order + 1; j += 2){
		for (int i=0; i < len; i++) 
			out_model_matrix[j + i*n_params] = cos(effdf*(j-poly_order+1)/2*i);
		for (int i=0; i < len; i++) 
			out_model_matrix[(j+1) + i*n_params] = sin(effdf*(j-poly_order+1)/2*i);
	}
	return m;
}



GslMatrixPtr make_notch_filter( double sample_frequency, const std::vector<double> & freqs, 
				int len){
	g3_assert(freqs.size() > 0);
	g3_assert(sample_frequency > 0);

	int n_params = 2 * freqs.size();
	GslMatrixPtr m = get_shared_gsl_matrix( len, n_params);
	double * out_model_matrix = m->data;

	for (size_t j=0; j < freqs.size(); j++){
		for (int i=0; i < len; i++){
			double trig_val = 2.0 * M_PI / sample_frequency * i * freqs[j];
			out_model_matrix[n_params * i + 2*j] = sin(trig_val);
			out_model_matrix[n_params * i + 2*j + 1] = cos(trig_val);
		}
	}
	return m;
}




GslMatrixPtr make_poly_abscissa_filter(G3TimestreamConstPtr abscissa, 
				       int poly_order, int index_poly_order,
				       int len){
	g3_assert(poly_order >= 0);
	if (index_poly_order < 0) index_poly_order = 0;

	G3Timestream ab( *abscissa );

	int n_params = 1 + poly_order + index_poly_order;
	GslMatrixPtr m = get_shared_gsl_matrix( len, n_params);
	double * out_model_matrix = m->data;

	double min_ab = ab[0];
	double max_ab = ab[0];
	for (size_t i = 1; i < ab.size(); i++) {
		if (ab[i] < min_ab) min_ab = ab[i];
		if (ab[i] > max_ab) max_ab = ab[i];
	}
	double delta = max_ab - min_ab;
	g3_assert(delta > 0);

	for (size_t i = 0; i < ab.size(); i++) {
		ab[i] = 2.0 * ( ab[i] - min_ab ) / delta - 1.0;
	}
	
	//poly or hpf always have a mean filter
	for (int i=0; i < len; i++){
		out_model_matrix[n_params * i] = 1.0f;
	} 

	//cover the first case
	if (poly_order >0){
		int j = 1;
		for (int i=0; i < len; i++){
			double x = ab[i];
			out_model_matrix[n_params * i + j] = x;
		}
	}
	
	for (int j = 2; j <  1 + poly_order; j++){
		for (int i=0; i < len; i++){
			double x = ab[i];
			double j_leg = j;
			out_model_matrix[n_params * i + j] = 
				((2*j_leg-1) * x * out_model_matrix[(j-1)+i*n_params] -
				 (j_leg-1) * out_model_matrix[(j-2)+i*n_params])/j_leg;
		}
	}


	if (index_poly_order > 0){
		int j = 1 + poly_order;
		for (int i=0; i < len; i++){
			double x = -1.0f + 2 * ((double)i)/((double) len - 1.0f);
			out_model_matrix[n_params * i + j] = x;
		}
	}
	
	for (int j = 2 + poly_order; j <  1 + index_poly_order + poly_order; j++){
		size_t off_index;
		if (j == 2 + poly_order) off_index = 0;
		else off_index = j - 2;

		for (int i=0; i < len; i++){
			double x = -1.0f + 2 * ((double)i)/((double) len - 1.0f);
			double j_leg = j - poly_order;
			out_model_matrix[n_params * i + j] = 
				((2*j_leg-1) * x * out_model_matrix[(j-1)+i*n_params] -
				 (j_leg-1) * out_model_matrix[(off_index) + i*n_params])/j_leg;
		}
	}


	return m;
}

#ifdef ACTUALLY_USE_LAPACK
void bulk_lss_fitter_lapacke(GslMatrixPtr model_matrix_in, GslMatrixPtr not_masked_model_matrix, 
			     const std::vector<double> & data,
			     std::vector<double> & output_data){

	std::vector<double>nd(data);

	lapack_int M = model_matrix_in->size1;
	lapack_int N = model_matrix_in->size2;

	gsl_matrix * model_matrix = gsl_matrix_alloc(M,N);
	gsl_matrix_memcpy(model_matrix, model_matrix_in.get());

	g3_assert(data.size() == M);

	lapack_int nrhs = 1;
	lapack_int lda = N;
	lapack_int ldb = nrhs;

	gsl_vector * y_fit = gsl_vector_alloc(M);
	gsl_vector * x = gsl_vector_alloc(N);

	LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',M,N,nrhs,model_matrix->data,lda,&(nd[0]),ldb);
	for (size_t i=0; i < N; i++){
		x->data[i] = nd[i];
	}
	//the beginning of nd now has the fit params
	gsl_blas_dgemv(CblasNoTrans, 1.0, not_masked_model_matrix.get(), x, 0, y_fit);
	for (size_t i=0; i < M; i++){
		output_data[i] = data[i] - y_fit->data[i];
	}
	gsl_vector_free(y_fit);
	gsl_vector_free(x);
	gsl_matrix_free(model_matrix);
}
#endif



void bulk_lss_fitter_qr(GslMatrixPtr model_matrix, GslMatrixPtr not_masked_model_matrix, 
			GslMatrixConstPtr qr, GslVectorPtr tau,
			const std::vector<double> & data,
			int len, int n_params, 
			std::vector<double> & output_data){
	g3_assert(len > n_params);
	output_data = data;
	
	log_debug("in bulk lls fitter part");
	gsl_vector gsl_v;
	gsl_v.size = output_data.size();
	gsl_v.stride = 1;
	gsl_v.data = &(output_data[0]);
	gsl_v.block = NULL;
	gsl_v.owner = 0;
	
	size_t M = qr->size1;
	size_t N = qr->size2;
	
	gsl_vector * x = gsl_vector_alloc(N);
	gsl_vector * resi = gsl_vector_alloc(M);
	gsl_vector * y_fit = gsl_vector_alloc(M);
	int status = gsl_linalg_QR_lssolve (qr.get(), tau.get(), &gsl_v, x, resi);
	gsl_blas_dgemv (CblasNoTrans, 1.0, not_masked_model_matrix.get(), x, 0, y_fit);
	for (size_t i=0; i < M; i++){
		gsl_v.data[i] -= y_fit->data[i];
	}
	gsl_vector_free(x);
	gsl_vector_free(resi);
	gsl_vector_free(y_fit);
}




void bulk_lss_fitter_svd(GslMatrixPtr model_matrix, GslMatrixPtr not_masked_model_matrix, 
			 GslMatrixConstPtr U, GslMatrixConstPtr V, GslVectorConstPtr S,
			 const std::vector<double> & data,
			 int len, int n_params, 
			 std::vector<double> & output_data){
	g3_assert(len > n_params);
	output_data = data;
	
	log_debug("in bulk lls fitter part");
	gsl_vector gsl_v;
	gsl_v.size = output_data.size();
	gsl_v.stride = 1;
	gsl_v.data = &(output_data[0]);
	gsl_v.block = NULL;
	gsl_v.owner = 0;
	
	size_t M = U->size1;
	size_t N = U->size2;
	
	gsl_vector * x = gsl_vector_alloc(N);
	gsl_vector * y_fit = gsl_vector_alloc(M);
	int status = gsl_linalg_SV_solve (U.get(), V.get(), S.get(), &gsl_v, x);
	gsl_blas_dgemv (CblasNoTrans, 1.0, not_masked_model_matrix.get(), x, 0, y_fit);
	for (size_t i=0; i < M; i++){
		gsl_v.data[i] -= y_fit->data[i];
	}
	gsl_vector_free(x);
	gsl_vector_free(y_fit);
}


void bulk_quick_cheap_fit(GslMatrixPtr model_matrix, GslMatrixPtr not_masked_model_matrix,
			  const std::vector<double> & data,
			  int len, int n_params,
			  std::vector<double> & output_data){
	g3_assert(len > n_params);
	output_data = data;

	size_t N = model_matrix->size2;

	log_debug("in bulk quick fitter part");
	gsl_vector gsl_v;
	gsl_v.size = output_data.size();
	gsl_v.stride = 1;
	gsl_v.data = &(output_data[0]);
	gsl_v.block = NULL;
	gsl_v.owner = 0;

	std::vector<double> selfps;
	std::vector<double> fs;
	std::vector<double> ls;

	for (int j=0; j < 2; j++){
		for (size_t i=0; i < N; i++){
			gsl_vector_view basis_fn = gsl_matrix_column(model_matrix.get(), i);
			gsl_vector_view unmasked_basis_fn = gsl_matrix_column(
							     not_masked_model_matrix.get(), i);
			if (j==0){
				double selfp;
				gsl_blas_ddot(&basis_fn.vector, &basis_fn.vector, &selfp);
				double f = gsl_vector_get(&basis_fn.vector, 0);
				double l = gsl_vector_get(&basis_fn.vector, len-1);
				selfp = ((selfp - (f*f)) + (selfp - (l*l))) / 2;

				selfps.push_back(selfp);
				fs.push_back(f);
				ls.push_back(l);
			}

			double dotp;
			gsl_blas_ddot(&basis_fn.vector, &gsl_v, &dotp);
			double first = fs[i] * gsl_vector_get(&gsl_v, 0);
			double last = ls[i] * gsl_vector_get(&gsl_v, len-1);
			dotp = -((dotp - first) + (dotp - last)) / (2 * selfps[i]);

			gsl_blas_daxpy(dotp, &unmasked_basis_fn.vector, &gsl_v);
		}
	}
}


void lls_filter_ts_data(const G3TimestreamMap & timestreams, 
			const GslMatrixPtr & model_matrix, 
			const FilterMask & mask, 
			G3TimestreamMap & out_timestreams){
	#ifdef USE_QR
		GslMatrixPtr unmasked_QR; 
		GslVectorPtr unmasked_Tau;
		get_qr_of_matrix(model_matrix,unmasked_QR,unmasked_Tau);
      	#elif USE_SVD
		GslMatrixPtr unmasked_U; 
		GslMatrixPtr unmasked_V; 
		GslVectorPtr unmasked_S;
		get_thin_svd_of_matrix(model_matrix,unmasked_U,unmasked_V,unmasked_S);
	#endif

	const unsigned int ts_len = timestreams.begin()->second->size();
	const unsigned int n_params = model_matrix->size2;
	
	std::vector<std::string> ts_lst(timestreams.size());
	size_t i=0;
	for (auto iter = timestreams.begin(); iter != timestreams.end(); iter++, i++){
		out_timestreams[iter->first]=G3TimestreamPtr(new G3Timestream(*(iter->second)));
		ts_lst[i] = iter->first;
	}
	
#ifdef LLS_OMP
#pragma omp parallel for
#endif
	for (i=0; i < ts_lst.size(); i++){
		GslMatrixPtr ts_model_mat;
		
		#ifdef USE_QR
		GslMatrixPtr QR;
		GslVectorPtr Tau;
		#elif USE_SVD
		GslMatrixPtr U;
		GslMatrixPtr V;
		GslVectorPtr S;
		#endif

		log_debug("checking mask");
		if ( mask.has_masked_pixels.find(ts_lst[i]) == mask.has_masked_pixels.end()){
			log_fatal("Attempting to poly filter timestream that is not"
				  " listed in the mask, ts id is: %s", ts_lst[i].c_str());
		}


		if (mask.has_masked_pixels.at(ts_lst[i])){
			//ts_model_mat is allocated in this
			make_mask_lls_model_matrix( model_matrix, mask.pixel_mask.at(ts_lst[i]),
						    n_params, ts_model_mat);
			#ifdef USE_QR
			get_qr_of_matrix(ts_model_mat,QR,Tau);
			#elif USE_SVD
			get_thin_svd_of_matrix(ts_model_mat,U,V,S);
			#endif
		}else{
			ts_model_mat = model_matrix;
			#ifdef USE_QR
			QR = unmasked_QR;
			Tau = unmasked_Tau;
			#elif USE_SVD
			U = unmasked_U;
			V = unmasked_V;
			S = unmasked_S;
			#endif
		}


		log_debug("subtracting");
#ifdef USE_QR
	bulk_lss_fitter_qr(ts_model_mat, model_matrix,	QR, Tau, 
			   *timestreams.at(ts_lst[i]), ts_len, n_params, 
			   *out_timestreams[ts_lst[i]]);

#elif USE_SVD
		bulk_lss_fitter_svd(ts_model_mat, model_matrix,	U,V,S, 
				    *timestreams.at(ts_lst[i]), ts_len, n_params, 
				    *out_timestreams[ts_lst[i]]);
#else
			bulk_quick_cheap_fit(ts_model_mat, model_matrix,
				     *timestreams.at(ts_lst[i]), ts_len, n_params,
				     *out_timestreams[ts_lst[i]]);
#endif
	}
}



G3TimestreamMapPtr poly_and_mhpf_filter_ts_data(G3TimestreamMapConstPtr timestreams, 
						FilterMaskConstPtr mask, 
						double sample_frequency, double freq_cutoff,
						int poly_order){
	G3TimestreamMapPtr out_timestreams(new G3TimestreamMap);
	int len = timestreams->at(timestreams->begin()->first)->size();
	GslMatrixPtr model_mat = make_poly_and_highpass_fit( sample_frequency,
							     freq_cutoff, poly_order, 
							     len);
	lls_filter_ts_data(*timestreams, model_mat, *mask, *out_timestreams);
	return out_timestreams;
}



G3TimestreamMapPtr poly_filter_ts_data_with_abscissa(
	G3TimestreamConstPtr abscissa,
	G3TimestreamMapConstPtr timestreams, 
	FilterMaskConstPtr mask, 
	int poly_order,
	int index_poly_order){

	G3TimestreamMapPtr out_timestreams(new G3TimestreamMap);
	int len = timestreams->at(timestreams->begin()->first)->size();
	GslMatrixPtr model_mat = 
		make_poly_abscissa_filter(abscissa, poly_order, index_poly_order, len);
	lls_filter_ts_data(*timestreams, model_mat, *mask, *out_timestreams);
	return out_timestreams;
}


G3TimestreamMapPtr notch_filter_lls(
	G3TimestreamMapConstPtr timestreams, 
	FilterMaskConstPtr mask, 
	const std::vector<double> & freqs){
	G3TimestreamMapPtr out_timestreams(new G3TimestreamMap);
	int len = timestreams->at(timestreams->begin()->first)->size();
	double sample_rate = timestreams->at(timestreams->begin()->first)->GetSampleRate();
	GslMatrixPtr model_mat = 
		make_notch_filter(sample_rate, freqs, len);
	lls_filter_ts_data(*timestreams, model_mat, *mask, *out_timestreams);
	return out_timestreams;
}


void rolling_mean_filter(const double * ts, int ts_len, int filter_width, double * out_ts){
	g3_assert(filter_width > 0);
	g3_assert(filter_width < ts_len);
	int filter_ind_low = filter_width/2;
	int filter_ind_high = filter_width/2 + filter_width%2;
	
	double sum_val = 0;
	double  ncount = 0;

	// sum values from 0 to filter_ind_high to get the 0th value
	for (int i=0; i < filter_ind_high; i++){
		if (std::isfinite(ts[i])){
			sum_val += ts[i];
			ncount++; 
		}
	}

	// start filtering and adding values to fill out our filter
	for (int i = 0; i < filter_ind_low; i++){
		out_ts[i] = ts[i] - sum_val/ncount; 
		if (std::isfinite(ts[i + filter_ind_high])){
			sum_val += ts[i + filter_ind_high];
			ncount++;
		}
	}

	for (int i = filter_ind_low; i < ts_len-filter_ind_high; i++){
		out_ts[i] = ts[i] - sum_val/ncount;
		if (std::isfinite(ts[i - filter_ind_low])){
			sum_val -= ts[i - filter_ind_low];
			ncount--;
		}
		if (std::isfinite(ts[i + filter_ind_high])){
			sum_val += ts[i + filter_ind_high];
			ncount++;
		}
	}
	for (int i = ts_len-filter_ind_high; i < ts_len; i++){
		out_ts[i] = ts[i] - sum_val/ncount;
		if (std::isfinite(ts[i - filter_ind_low])){
			sum_val -= ts[i - filter_ind_low];
			ncount--;
		}
	}
}




void mean_smooth_filter(const double * ts, int ts_len, int filter_width, double * out_ts){
	g3_assert(filter_width > 0);
	g3_assert(filter_width < ts_len);

	int filter_ind_low = filter_width/2;
	int filter_ind_high = filter_width/2 + filter_width%2;
	
	double sum_val = 0;
	double  ncount = 0;

	// sum values from 0 to filter_ind_high to get the 0th value
	for (int i=0; i < filter_ind_high; i++){
		sum_val += ts[i];
		ncount++; 
	}

	// start filtering and adding values to fill out our filter
	for (int i = 0; i < filter_ind_low; i++){
		out_ts[i] = sum_val/ncount; 
		sum_val += ts[i + filter_ind_high];
		ncount++;
	}

	for (int i = filter_ind_low; i < ts_len-filter_ind_high; i++){
		sum_val += -1* ts[i-filter_ind_low] + ts[i+filter_ind_high];
		out_ts[i] = sum_val/ncount;
	}

	for (int i = ts_len-filter_ind_high; i < ts_len; i++){
		sum_val -= ts[i - filter_ind_low];
		ncount--;
		out_ts[i] = sum_val/ncount;
	}
}


G3TimestreamPtr mean_smooth_ts(G3TimestreamConstPtr  timestreams, int filter_width){
	G3TimestreamPtr out_timestreams = G3TimestreamPtr(new G3Timestream(*timestreams));
	mean_smooth_filter((&(*timestreams)[0]), out_timestreams->size(), 
			    filter_width, &((*out_timestreams)[0]));
	return out_timestreams;
}




G3TimestreamPtr rolling_mean_filter_ts(G3TimestreamConstPtr  timestreams, int filter_width){
	G3TimestreamPtr out_timestreams = G3TimestreamPtr(new G3Timestream(*timestreams));
	rolling_mean_filter((&(*timestreams)[0]), out_timestreams->size(), 
			    filter_width, &((*out_timestreams)[0]));
	return out_timestreams;
}

void rolling_mean_filter_ts_data(G3TimestreamMapConstPtr  timestreams, int filter_width,
				 G3TimestreamMapPtr  out_timestreams){
	out_timestreams = G3TimestreamMapPtr(new G3TimestreamMap);
	for (auto it = timestreams->begin(); it != timestreams->end(); it++){
		(*out_timestreams)[it->first] = G3TimestreamPtr(new G3Timestream(*it->second));
		rolling_mean_filter((&(*it->second)[0]), out_timestreams->at(it->first)->size(), 
				    filter_width, &(*out_timestreams->at(it->first))[0]);
	}
}


namespace bp = boost::python;
PYBINDINGS("todfilter"){
	EXPORT_FRAMEOBJECT(FilterMask, init<>(), 
			   "This stores the information about whether an individual"
			   " sample for an individual detector is masked.  1 (true)"
			   " means a sample for a detector is masked (masked means"
			   " not used in the fit), 0 (false) means it is not.  It"
			   " contains two maps, pixel_mask actually specifies the masked"
			   " samples.  has_masked_pixels specifies whether any values"
			   " in the pixel_mask are true for a specific key."
		)
		.def_readwrite("has_masked_pixels", 
			       &FilterMask::has_masked_pixels, 
			       "A G3MapInt with the same keys as pixel_mask.  "
			       "The value mapped to should be equal to applying 'or'"
			       " to all the values in pixel_mask.  It says whether"
			       " any pixels are masked in the array.")
		.def_readwrite("pixel_mask", &FilterMask::pixel_mask, "A G3MapVectorInt"
			       " that maps whether or not a sample is masked when doing"
			       " linear least squares fitting.  1 means the pixel is"
			       " masked (not used in LLS model fitting). The keys are bolometer ids.")
		.def("is_valid", &FilterMask::is_valid)
		;

	register_pointer_conversions<FilterMask>();

	bp::def("make_empty_filter_mask", make_empty_filter_mask, bp::arg("ts_map"),
		"Returns an empty FilterMask for the timestream map."
		);

	bp::def("make_filter_mask", make_filter_mask,
		(bp::arg("point_source_mask,"), bp::arg("detector_pointing")),
		"Generates a filter mask from the detector pointing information and the "
		"point source mask map.  It returns the filter mask."
		);

	bp::def("poly_and_mhpf_filter_ts_data", poly_and_mhpf_filter_ts_data,
		(bp::arg("timestreams"), bp::arg("mask"), 
		 bp::arg("sample_frequency"), bp::arg("freq_cutoff"), 
		 bp::arg("poly_order")),
		"Returns a G3TimestreamMap of the filtered timestreams.  Does masked linear "
		"least squares fitting of the timestream"
		);
    
    bp::def("notch_filter_lls", notch_filter_lls,
        (bp::arg("timestreams"), bp::arg("mask"), 
         bp::arg("freqs")),
        "Returns a G3TimestreamMap of the filtered timestreams. "
        "The notching is based on masked linear least squares fitting of the timestream.");

	bp::def("get_polarization_fit_params", get_polarization_fit_params,
		(bp::arg("Q"), bp::arg("U"), bp::arg("mask"), bp::arg("map")),
		"This code uses linear least squares to fit the Q and U maps.  It returns "
		"the two fit parameters, Q, U.  The mask determines the area used in the fit."
		"In the mask, 1 means it's used, 0 means it is not used.\n\n"
		"The map used should be (individual detector map - coadd T map)."
		);

	bp::def("rolling_mean_filter_ts_data", rolling_mean_filter_ts_data,
		(bp::arg("timestreams"), bp::arg("filter_width"), bp::arg("out_timestreams"))
		);

	bp::def("rolling_mean_filter_ts", rolling_mean_filter_ts, 
		( bp::arg("timestream"), bp::arg("filter_width")),
		"Returns the rolling mean filtered timestream."
		);

	bp::def("mean_smooth_ts", mean_smooth_ts, 
		(bp::arg("timestreams"), bp::arg("filter_width")),
		"Convolves the timestreams with a tophat of width filter_width.");

	bp::def("poly_filter_ts_data_with_abscissa", poly_filter_ts_data_with_abscissa,
		(bp::arg("abscissa"), bp::arg("timestreams"), bp::arg("mask"), 
		 bp::arg("poly_order"), bp::arg("index_poly_order")
			),
		"Does a masked linear least squares filter of timestreams, where you can"
		"An x value (abscissa) to the polynomial filter.  While fitting the data "
		"It uses this and at the same time does a polynomial filter in time.  ",
		"poly_order is the order of the fit with the abscissa.  index_poly_order is the"
		"order the polynomial order for the time based fit."
		);

	bp::def("simple_common_mode_filter_cpp", simple_common_mode_filter_cpp);
	bp::def("masked_common_mode_filter_cpp", masked_common_mode_filter_cpp);
	bp::def("notch_filter_lls", notch_filter_lls,
		(bp::arg("timestreams"), bp::arg("mask"), bp::arg("notch_frequencies")),
		"Fits and then subtracts off sine and cosine waves with frequencies "\
		"in the vector notch_frequencies to the timestreams in timestreams. "\
		"mask is the FilterMask.\n\nReturns a new filtered timestream map\n\n"\
		"Because this is fitting one frequency of sine wave you may need to "\
		"filter multiple sine waves close in frequency space to remove a finite "\
		"bandwidth."
		);
}

