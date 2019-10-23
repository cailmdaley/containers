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
#define max(a,b) ((a < b) ? (b) : (a))

static PyObject *tmplfilter(PyObject *self, PyObject *args);

static PyMethodDef methods[] = {
        { "tmplfilter", tmplfilter, METH_VARARGS },
        { NULL, NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef filter_module = {
	PyModuleDef_HEAD_INIT, "tmplfilter", "Process map with a real-space filter", -1, methods
};

PyObject *PyInit_tmplfilter(void)
{
	PyObject *module = PyModule_Create(&filter_module);

#else
void inittmplfilter(void)
{
#endif

	import_array();
	
#if PY_MAJOR_VERSION >= 3
	return module;
#else
	Py_InitModule("tmplfilter", methods);
#endif
}

static PyObject *tmplfilter(PyObject *self, PyObject *args)
{
	double *map, *errors, *filter;
	double *outmap, *outerrors;
	PyObject *xmap, *xerrors, *xfilter;
	PyArrayObject *amap, *aerrors, *afilter;
	PyArrayObject *aoutmap, *aouterrors;
	int i, j, k, m;
	double num, dem;
	double mval, fval, eval;

	if (!PyArg_ParseTuple(args, "OOO", &xmap, &xerrors, &xfilter))
		return NULL;

	amap = (PyArrayObject *)PyArray_ContiguousFromObject(xmap,
	    NPY_DOUBLE, 2, 2);
	aerrors = (PyArrayObject *)PyArray_ContiguousFromObject(xerrors,
	    NPY_DOUBLE, 2, 2);
	afilter = (PyArrayObject *)PyArray_ContiguousFromObject(xfilter, NPY_DOUBLE,
	    2, 2);
	map = PyArray_DATA(amap);
	errors = PyArray_DATA(aerrors);
	filter = PyArray_DATA(afilter);

	if (PyArray_DIM(afilter, 0) % 2 == 0 || PyArray_DIM(afilter, 1) % 2 == 0) {
		PyErr_SetString(PyExc_ValueError, "Filter needs to have an odd number of pixels");
		return NULL;
	}

	/*
	 * The following fits a copy of the filtering kernel centered at every map
	 * pixel to the input map using linear least-squares. The output is the
	 * multiplicative scaling factor of the best-fit filtering kernel and the
	 * error from this fit.
	 */

	aoutmap = (PyArrayObject *)PyArray_SimpleNew(2, PyArray_DIMS(amap), NPY_DOUBLE);
	aouterrors = (PyArrayObject *)PyArray_SimpleNew(2, PyArray_DIMS(amap), NPY_DOUBLE);
	outmap = PyArray_DATA(aoutmap);
	outerrors = PyArray_DATA(aouterrors);
	for (i = 0; i < PyArray_DIM(amap, 0); i++) {
		for (j = 0; j < PyArray_DIM(amap, 1); j++) {
			num = dem = 0;
			for (k = -PyArray_DIM(afilter, 0)/2; k <= PyArray_DIM(afilter, 0)/2; k++) {
				if (i + k < 0 || i + k >= PyArray_DIM(amap, 0))
					continue;
				for (m = -PyArray_DIM(afilter, 1)/2; m <= PyArray_DIM(afilter, 1)/2; m++) {
					if (j + m < 0 || j + m >= PyArray_DIM(amap, 1))
						continue;

					eval = errors[(i+k)*PyArray_STRIDES(amap)[0]/sizeof(double) + j+m];
					fval = filter[(k + PyArray_DIM(afilter, 0)/2)*PyArray_STRIDES(afilter)[0]/sizeof(double) + m + PyArray_DIM(afilter, 1)/2];
					if (!isfinite(eval) || eval == 0 || fval == 0)
						continue;
					mval = map[(i+k)*PyArray_STRIDES(amap)[0]/sizeof(double) + j+m];

					num += mval*fval/pow(eval, 2);
					dem += pow(fval, 2)/pow(eval, 2);
				}
			}
			outmap[i*PyArray_STRIDES(aoutmap)[0]/sizeof(double) + j] = num/dem;
			outerrors[i*PyArray_STRIDES(aoutmap)[0]/sizeof(double) + j] = 1./sqrt(dem);
		}
	}

	Py_DECREF(amap);
	Py_DECREF(aerrors);
	Py_DECREF(afilter);

	return Py_BuildValue("OO", aoutmap, aouterrors);
}

