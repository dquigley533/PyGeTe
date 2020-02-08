/* -*- mode: C */

%define DOCSTRING
"Implements total energy and atom local energy calculations for
atomistic simulations of germanium telluride using the model of
Zipoli and Curioni, New J. Phys. 15, 123006 (2013)." 
%enddef

/* energy.i */
/* N.B. Implementing module docstring using the method described 
at http://swig.org/Doc1.3/Python.html#Python_nn66 stops distutils
from recognising the module name.... */
%module energy
%{
#define SWIG_FILE_WITH_INIT

/* This will be included in the wrapper code */
#include "energy.h"
%}

/* Numpy array typemap */
%include "numpy.i"

%init %{
  import_array();
%}

/* Map array array onto arguments used in the C interface */

/* Atom positions */
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int n, int d, double* pos)};

/* Matrix of cell vectors */
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int dh2, int dh1, double* hmatrix)};

/* Species array */
%apply (int DIM1, int* IN_ARRAY1) {(int n2, int* species)};

%feature("docstring", "This is a docstring for compute_neighbour_list()") compute_neighbour_list;


/* This will be parsed to generate the wrapper */
%include "energy.h"

