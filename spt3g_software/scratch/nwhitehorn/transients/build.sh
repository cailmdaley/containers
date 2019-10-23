python_includes=$SROOT/include/python?.?/
numpy_includes=$SROOT/lib/python?.?/site-packages/numpy/core/include
set -x
cc -fPIC -shared -o tmplfilter.so -I $SROOT/include -I $python_includes -I $numpy_includes tmplfilter.c
