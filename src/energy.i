/* -*- mode: C */

%define DOCSTRING
"Implements total energy and atom local energy calculations for
atomistic simulations of germanium telluride using the model of
Zipoli and Curioni, New J. Phys. 15, 123006 (2013)." 
%enddef

/* energy.i */
/* N.B. Implementing module docstring using the method described 
at http://swig.org/Doc1.3/Python.html#Python_nn66 stops distutils
from recognising the module name.... 
%module(docstrig=DOCSTRING) energy
*/
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

/* Docstring information for compute_neighbour_list */
%feature("autodoc", "compute_neighbour_list(positions, cellmatrix, species)") compute_neighbour_list;
%define cpt_nl_string
"
    Updates the internal list of atoms that lie within (cutoff+skin) of each atom. 

    This function should be called whenever atoms may have moved by a sufficient 
    distance than any previous neighbour list will no longer be valid. The cutoff
    distance defined by the model potential is species pair dependent. The skin
    width is fixed at 1.0 Angstrom.

    Parameters
    ----------

    positions  : Numpy array of size (N, 3). Compatible with the positions attribute
                 of an ASE atoms object. Holds atomic position vectors.

    cellmatrix : Numpy array of size (3,3). Holds the 3 vectors defining the 
                 parallelepiped simulation cell. Compatible with the cell attribute
                 of an ASE atoms object.

    species    : List of N atoms species types defined using the constants GE=0, TE=1. 


"
%enddef
%feature("docstring", cpt_nl_string) compute_neighbour_list;




/* This will be parsed to generate the wrapper */
%include "energy.h"

