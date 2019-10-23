#include <Python.h>
#include <numpy/numpyconfig.h>
#ifdef NPY_1_7_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif
#include <numpy/ndarrayobject.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_psi.h>
#include <stdio.h>

#define min(a,b) ((a < b) ? (a) : (b))

static PyObject *pyfitlightcurve(PyObject *self, PyObject *args);
static PyObject *pyllheval(PyObject *self, PyObject *args);
static double fitlightcurve(int n, double *points, double *errors,
    int ntemplate, double *template, double *times, double *amp, double *t0,
    double *sigma, double tstep);

static PyMethodDef methods[] = {
        { "fitlightcurve", pyfitlightcurve, METH_VARARGS },
        { "llheval", pyllheval, METH_VARARGS },
        { NULL, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef fitlc_module = {
	PyModuleDef_HEAD_INIT, "fitlc", "fitlc module", -1, methods
};

PyObject *PyInit_fitlc(void)
{
	PyObject *module = PyModule_Create(&fitlc_module);

#else
void initfitlc(void)
{
#endif

	import_array();
	
#if PY_MAJOR_VERSION >= 3
	return module;
#else
	Py_InitModule("fitlc", methods);
#endif
}

static PyObject *pyfitlightcurve(PyObject *self, PyObject *args)
{
	int n, ntempl, i;
	double bestllhdiff, bestt0;
	double *points, *errors, *t, *templ;
	double llhdiff, amp, t0, sigma;
	PyObject *xpoints, *xerrors, *xt, *xtempl;
	PyArrayObject *apoints, *aerrors, *at, *atempl;

	if (!PyArg_ParseTuple(args, "OOOO", &xpoints, &xerrors, &xt, &xtempl))
		return NULL;

	apoints = (PyArrayObject *)PyArray_ContiguousFromObject(xpoints,
	    NPY_DOUBLE, 1, 1);
	aerrors = (PyArrayObject *)PyArray_ContiguousFromObject(xerrors,
	    NPY_DOUBLE, 1, 1);
	at = (PyArrayObject *)PyArray_ContiguousFromObject(xt, NPY_DOUBLE,
	    1, 1);
	atempl = (PyArrayObject *)PyArray_ContiguousFromObject(xtempl,
	    NPY_DOUBLE, 1, 1);
	n = PyArray_DIM(apoints, 0);
	points = PyArray_DATA(apoints);
	errors = PyArray_DATA(aerrors);
	t = PyArray_DATA(at);
	ntempl = PyArray_DIM(atempl, 0);
	templ = PyArray_DATA(atempl);
	if (n != PyArray_DIM(aerrors, 0)) {
		PyErr_SetString(PyExc_ValueError, "Data and error array "
		    "dimensions do not match");
		return NULL;
	}
	if (n != PyArray_DIM(at, 0)) {
		PyErr_SetString(PyExc_ValueError, "Data and time array "
		    "dimensions do not match");
		return NULL;
	}
	if (ntempl < 1) {
		PyErr_SetString(PyExc_ValueError, "Too few template points");
		return NULL;
	}
	if (ntempl > n) {
		PyErr_SetString(PyExc_ValueError, "Too many template points");
		return NULL;
	}
	if ((n % ntempl) != 0) {
		PyErr_SetString(PyExc_ValueError, "Data points not an even "
		    "multiple of template points");
		return NULL;
	}
	
	bestt0 = 0;
	bestllhdiff = 0;
	#define NSEEDS 30.
	for (i = 0; i <= NSEEDS; i++) {
		/* Try ten different seeds to check for convergence */
		amp = 0; t0 = t[0] + (t[n-1] - t[0])*i/NSEEDS; sigma = 10;
		llhdiff = fitlightcurve(n, points, errors, ntempl, templ, t,
		    &amp, &t0, &sigma, 1.5*(t[n-1] - t[0])/NSEEDS);

		t0 = t[0] + (t[n-1] - t[0])*i/NSEEDS;
		if (llhdiff < bestllhdiff) {
			bestllhdiff = llhdiff;
			bestt0 = t0;
		}
	}
	amp = 0; t0 = bestt0; sigma = 10;
	llhdiff = fitlightcurve(n, points, errors, ntempl, templ, t,
	    &amp, &t0, &sigma, 1.5*(t[n-1] - t[0])/NSEEDS);

	Py_DECREF(apoints);
	Py_DECREF(aerrors);
	Py_DECREF(at);
	Py_DECREF(atempl);

	return Py_BuildValue("dddd", llhdiff, amp, t0, sigma);
}

struct dataparams {
	double *x;
	double *t;
	double *errors;
	int n;

	double *template;
	int ntemplate;

	double *templ1, *templ2, *templ3;
};

static void gaussresid(const gsl_vector *x, void *xparams, double *f, gsl_vector *g)
{
	struct dataparams *params = xparams;
	double norm, t0, sigma;
	double dnorm = 0, dt0 = 0, dsigma = 0;
	double t, val;
	double llh = 0;
	int i, j;

	t0 = gsl_vector_get(x, 0);
	sigma = fabs(gsl_vector_get(x, 1));
	norm = gsl_vector_get(x, 2);

	if (g)
		gsl_vector_set_all(g, 0);
	
	for (i = 0; i < params->n/params->ntemplate; i++) {
		t = params->t[i*params->ntemplate];

		val = exp(-pow((t - t0)/(M_SQRT2*sigma), 2));
		if (!isfinite(val))
			fprintf(stderr, "Non-finite likelihood (%e), t %e t0 %e sigma %e norm %e\n", val, t, t0, sigma, norm);
		assert(isfinite(val));

		llh += -2*val*norm*params->templ2[i] + val*val*norm*norm*params->templ3[i];

		if (g) {
			dnorm = -2*val*params->templ2[i] + val*val*2*norm*params->templ3[i];
			val *= norm;
			
			dt0 = val*(t - t0)/pow(sigma, 2);
			dsigma = val*pow(t - t0, 2)/pow(sigma, 3);

			gsl_vector_set(g, 0, gsl_vector_get(g, 0) + -2*dt0*params->templ2[i] + 2*val*dt0*params->templ3[i]);
			gsl_vector_set(g, 1, gsl_vector_get(g, 1) + -2*dsigma*params->templ2[i] + 2*val*dsigma*params->templ3[i]);
			gsl_vector_set(g, 2, gsl_vector_get(g, 2) + dnorm);
		}
	}

	assert(isfinite(llh));

	if (t0 < 0) {
		llh += t0*t0;
		if (g)
			gsl_vector_set(g, 0, gsl_vector_get(g, 0) + 2*t0);
	}
	assert(isfinite(llh));
	if (t0 > params->t[params->n-1]) {
		llh += pow(t0 - params->t[params->n-1], 2);
		if (g)
			gsl_vector_set(g, 0, gsl_vector_get(g, 0) + 2*(t0 - params->t[params->n-1] + 10));
	}
	assert(isfinite(llh));

	/* Apply flare width penalty */
	if (1)
		llh -= 2*(log(min(2.355*fabs(sigma), params->t[params->n-1] - params->t[0])) - log(params->t[params->n-1] - params->t[0]));
	assert(isfinite(llh));

	if (f)
		*f = llh;
}

static double gaussresid_f(const gsl_vector *x, void *xparams)
{
	double ret;

	gaussresid(x, xparams, &ret, NULL);
	return (ret);
}

static void gaussresid_df(const gsl_vector *x, void *xparams, gsl_vector *g)
{

	gaussresid(x, xparams, NULL, g);
}

static PyObject *pyllheval(PyObject *self, PyObject *args)
{
	int n, ntempl, i, j;
	double bestllhdiff, bestt0;
	double *points, *errors, *t, *templ;
	double llhdiff, amp, t0, sigma;
	gsl_vector *params;
	PyObject *xpoints, *xerrors, *xt, *xtempl;
	PyArrayObject *apoints, *aerrors, *at, *atempl;
	struct dataparams data;

	if (!PyArg_ParseTuple(args, "OOOOddd", &xpoints, &xerrors, &xt, &xtempl, &amp, &t0, &sigma))
		return NULL;

	apoints = (PyArrayObject *)PyArray_ContiguousFromObject(xpoints,
	    NPY_DOUBLE, 1, 1);
	aerrors = (PyArrayObject *)PyArray_ContiguousFromObject(xerrors,
	    NPY_DOUBLE, 1, 1);
	at = (PyArrayObject *)PyArray_ContiguousFromObject(xt, NPY_DOUBLE,
	    1, 1);
	atempl = (PyArrayObject *)PyArray_ContiguousFromObject(xtempl,
	    NPY_DOUBLE, 1, 1);
	n = PyArray_DIM(apoints, 0);
	points = PyArray_DATA(apoints);
	errors = PyArray_DATA(aerrors);
	t = PyArray_DATA(at);
	ntempl = PyArray_DIM(atempl, 0);
	templ = PyArray_DATA(atempl);
	if (n != PyArray_DIM(aerrors, 0)) {
		PyErr_SetString(PyExc_ValueError, "Data and error array "
		    "dimensions do not match");
		return NULL;
	}
	if (n != PyArray_DIM(at, 0)) {
		PyErr_SetString(PyExc_ValueError, "Data and time array "
		    "dimensions do not match");
		return NULL;
	}
	if (ntempl < 1) {
		PyErr_SetString(PyExc_ValueError, "Too few template points");
		return NULL;
	}
	if (ntempl > n) {
		PyErr_SetString(PyExc_ValueError, "Too many template points");
		return NULL;
	}
	if ((n % ntempl) != 0) {
		PyErr_SetString(PyExc_ValueError, "Data points not an even "
		    "multiple of template points");
		return NULL;
	}

	data.x = points;
	data.errors = errors;
	data.t = t;
	data.n = n;
	data.template = templ;
	data.ntemplate = ntempl;

	data.templ1 = calloc(n/ntempl, sizeof(double));
	data.templ2 = calloc(n/ntempl, sizeof(double));
	data.templ3 = calloc(n/ntempl, sizeof(double));

	for (i = 0; i < n; i+= ntempl) {
		for (j = 0; j < ntempl; j++) {
			if (!isfinite(points[i+j]) || errors[i+j] == 0 ||
			    !isfinite(errors[i+j]))
				continue;

			data.templ1[i/ntempl] += pow(points[i+j],2)/pow(errors[i+j],2);
			data.templ2[i/ntempl] += points[i+j]*templ[j]/pow(errors[i+j],2);
			data.templ3[i/ntempl] += pow(templ[j],2)/pow(errors[i+j],2);
		}
	}

	params = gsl_vector_alloc(3);
	gsl_vector_set(params, 0, t0);
	gsl_vector_set(params, 1, sigma);
	gsl_vector_set(params, 2, amp);
	llhdiff = gaussresid_f(params, &data);
	gsl_vector_free(params);

	free(data.templ1);
	free(data.templ2);
	free(data.templ3);

	Py_DECREF(apoints);
	Py_DECREF(aerrors);
	Py_DECREF(at);
	Py_DECREF(atempl);

	return Py_BuildValue("d", llhdiff);
}

static double fitlightcurve(int n, double *points, double *errors,
    int ntemplate, double *template, double *times, double *amp, double *t0,
    double *sigma, double tstep)
{
	gsl_multimin_fminimizer *solver;
	gsl_multimin_function f;
	gsl_vector *params, *ss;
	struct dataparams data;
	int i, j, iter, status;
	double llhdiff;

	data.x = points;
	data.errors = errors;
	data.t = times;
	data.n = n;
	data.template = template;
	data.ntemplate = ntemplate;

	data.templ1 = calloc(n/ntemplate, sizeof(double));
	data.templ2 = calloc(n/ntemplate, sizeof(double));
	data.templ3 = calloc(n/ntemplate, sizeof(double));

	for (i = 0; i < n; i+= ntemplate) {
		for (j = 0; j < ntemplate; j++) {
			if (!isfinite(points[i+j]) || errors[i+j] == 0 ||
			    !isfinite(errors[i+j]))
				continue;

			data.templ1[i/ntemplate] += pow(points[i+j],2)/pow(errors[i+j],2);
			data.templ2[i/ntemplate] += points[i+j]*template[j]/pow(errors[i+j],2);
			data.templ3[i/ntemplate] += pow(template[j],2)/pow(errors[i+j],2);
		}
	}

	/* Solve for gauss fit */
	params = gsl_vector_alloc(3);
	ss = gsl_vector_alloc(3);
	solver = gsl_multimin_fminimizer_alloc(
	    gsl_multimin_fminimizer_nmsimplex2rand, 3);
	gsl_vector_set(params, 0, *t0);
	gsl_vector_set(params, 1, tstep/1.5 /* sigma */);
	gsl_vector_set(params, 2, 0 /* XXX calc? */);
	gsl_vector_set(ss, 0, tstep);
	gsl_vector_set(ss, 1, tstep);
	gsl_vector_set(ss, 2, 100);
	f.f = &gaussresid_f;
	f.n = 3;
	f.params = &data;
	gsl_multimin_fminimizer_set(solver, &f, params, ss);

	iter = 0;
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(solver);
		if (status)
			break;
		status = gsl_multimin_test_size(gsl_multimin_fminimizer_size(solver), 1e-2);
	} while (status == GSL_CONTINUE && iter < 1000);

	llhdiff = solver->fval;
	if (status != GSL_SUCCESS)
		llhdiff = INFINITY;

	*t0 = gsl_vector_get(solver->x, 0);
	*sigma = fabs(gsl_vector_get(solver->x, 1));
	*amp = gsl_vector_get(solver->x, 2);
	gsl_vector_free(params);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(solver);

	free(data.templ1);
	free(data.templ2);
	free(data.templ3);

	return (llhdiff);
}

