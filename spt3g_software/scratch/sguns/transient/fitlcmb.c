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
//static PyObject *pyllheval(PyObject *self, PyObject *args);
static double fitlightcurve(int n, double *points, double *errors,
    int nbands, double *times, double *amp, double *t0,
    double *sigma, double tstep);

static PyMethodDef methods[] = {
        { "fitlightcurve", pyfitlightcurve, METH_VARARGS },
//        { "llheval", pyllheval, METH_VARARGS },
        { NULL, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef fitlcmb_module = {
	PyModuleDef_HEAD_INIT, "fitlcmb", "fitlcmb module", -1, methods
};

PyObject *PyInit_fitlcmb(void)
{
	PyObject *module = PyModule_Create(&fitlcmb_module);

#else
void initfitlcmb(void)
{
#endif

	import_array();
	
#if PY_MAJOR_VERSION >= 3
	return module;
#else
	Py_InitModule("fitlcmb", methods);
#endif
}

static PyObject *pyfitlightcurve(PyObject *self, PyObject *args)
{
	int n, nbands, i, j;
	double bestllhdiff, bestt0;
	double *points, *errors, *t;
	double llhdiff, t0, sigma;
	PyObject *xpoints, *xerrors, *xt;
	PyArrayObject *apoints, *aerrors, *at;

	if (!PyArg_ParseTuple(args, "OOOi", &xpoints, &xerrors, &xt, &nbands))
		return NULL;

	apoints = (PyArrayObject *)PyArray_ContiguousFromObject(xpoints,
	    NPY_DOUBLE, 1, 1);
	aerrors = (PyArrayObject *)PyArray_ContiguousFromObject(xerrors,
	    NPY_DOUBLE, 1, 1);
	at = (PyArrayObject *)PyArray_ContiguousFromObject(xt, NPY_DOUBLE,
	    1, 1);
	n = PyArray_DIM(apoints, 0);
	points = PyArray_DATA(apoints);
	errors = PyArray_DATA(aerrors);
	t = PyArray_DATA(at);

        double *amp;
        amp = malloc(nbands * sizeof(double));

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
	
	bestt0 = 0;
	bestllhdiff = 0;
	#define NSEEDS 30.
	for (i = 0; i <= NSEEDS; i++) {
		/* Try ten different seeds to check for convergence */
                for(j = 0; j < nbands; j++) {
                    *(amp + j) = 0.;
                }
		t0 = t[0] + (t[n-1] - t[0])*i/NSEEDS; sigma = 10;
		llhdiff = fitlightcurve(n, points, errors, nbands, t,
		    amp, &t0, &sigma, 1.5*(t[n-1] - t[0])/NSEEDS);

		t0 = t[0] + (t[n-1] - t[0])*i/NSEEDS;
		if (llhdiff < bestllhdiff) {
			bestllhdiff = llhdiff;
			bestt0 = t0;
		}
	}
        for(j = 0; j < nbands; j++) {
            *(amp + j) = 0.;
        }
	t0 = bestt0; sigma = 10;
	llhdiff = fitlightcurve(n, points, errors, nbands, t,
	    amp, &t0, &sigma, 1.5*(t[n-1] - t[0])/NSEEDS);

	Py_DECREF(apoints);
	Py_DECREF(aerrors);
	Py_DECREF(at);

        PyObject *lamp = PyTuple_New(nbands);
        for (j = 0; j < nbands;j++) {
            PyObject *num = PyFloat_FromDouble(amp[j]);
            if (!num) {
                Py_DECREF(lamp);
                return NULL;
            }
            PyTuple_SET_ITEM(lamp, j, num);
        }

	return Py_BuildValue("dOdd", llhdiff, lamp, t0, sigma);
}

struct dataparams {
	double *x;
	double *t;
	double *errors;
	int n;
        int nbands;
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

	if (g)
		gsl_vector_set_all(g, 0);
	
	for (i = 0; i < params->n/params->nbands; i++) {
		t = params->t[params->nbands*i];

		val = exp(-pow((t - t0)/(M_SQRT2*sigma), 2));
		if (!isfinite(val))
			fprintf(stderr, "Non-finite likelihood (%e), t %e t0 %e sigma %e\n", val, t, t0, sigma);
		assert(isfinite(val));

                for (j = 0; j < params->nbands; j++) {
                    norm = gsl_vector_get(x, 2+j);
                    llh += -2*val*norm*params->x[params->nbands*i+j]/pow(params->errors[params->nbands*i+j],2) + val*val*norm*norm/pow(params->errors[params->nbands*i+j],2);
                }

		if (g) {
                        fprintf(stderr, "This mode is not supported for now");
                        assert(0);
//			dnorm = -2*val*params->templ2[i] + val*val*2*norm*params->templ3[i];
//			val *= norm;
//			
//			dt0 = val*(t - t0)/pow(sigma, 2);
//			dsigma = val*pow(t - t0, 2)/pow(sigma, 3);
//
//			gsl_vector_set(g, 0, gsl_vector_get(g, 0) + -2*dt0*params->templ2[i] + 2*val*dt0*params->templ3[i]);
//			gsl_vector_set(g, 1, gsl_vector_get(g, 1) + -2*dsigma*params->templ2[i] + 2*val*dsigma*params->templ3[i]);
//			gsl_vector_set(g, 2, gsl_vector_get(g, 2) + dnorm);
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

static double fitlightcurve(int n, double *points, double *errors,
    int nbands, double *times, double *amp, double *t0,
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
        data.nbands = nbands;

	params = gsl_vector_alloc(2+nbands);
	ss = gsl_vector_alloc(2+nbands);
	solver = gsl_multimin_fminimizer_alloc(
	    gsl_multimin_fminimizer_nmsimplex2rand, 2+nbands);
	gsl_vector_set(params, 0, *t0);
	gsl_vector_set(params, 1, tstep/1.5 /* sigma */);
	gsl_vector_set(ss, 0, tstep);
	gsl_vector_set(ss, 1, tstep);
        for (j = 0; j < nbands; j++) {
            gsl_vector_set(params, 2+j, 0 /* XXX calc? */);
            gsl_vector_set(ss, 2+j, 100);
        }
	f.f = &gaussresid_f;
	f.n = 2+nbands;
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
        for (j = 0; j < nbands; j++) {
	    *(amp+j) = gsl_vector_get(solver->x, 2+j);
        }
	gsl_vector_free(params);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(solver);

	return (llhdiff);
}

