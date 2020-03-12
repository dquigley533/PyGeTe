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

/* Docstring information for compute_neighbour_list */
%feature("autodoc", "compute_model_energy(positions, cellmatrix, species)") compute_model_energy;
%define cpt_mdl_string
"
    Compute the total energy of the system.

    Computes the total energy of the system (using the current neighbour list). To be
    used whenever the total energy of the system is required, e.g. when computing the
    change in energy due to a change in simulation cell shape and/or volume.

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
%feature("docstring", cpt_mdl_string) compute_model_energy;

%feature("autodoc", "compute_local_real_energy(iatom, positions, cellmatrix, species)") compute_local_real_energy;
%define cpt_loc_string
"
    Compute the contribution to total energy involving atom iatom.

    Computes all terms in the model potential (using the current neighbour list) 
    which depend on the position of a single atom. To be used when making single particle
    moves and hence only the energy contribution of the moved atom is required. 

    Parameters
    ----------

    iatom      : Integer index of the atom for which the energy contribution is required.

    positions  : Numpy array of size (N, 3). Compatible with the positions attribute
                 of an ASE atoms object. Holds atomic position vectors.

    cellmatrix : Numpy array of size (3,3). Holds the 3 vectors defining the 
                 parallelepiped simulation cell. Compatible with the cell attribute
                 of an ASE atoms object.

    species    : List of N atoms species types defined using the constants GE=0, TE=1. 


"
%enddef
%feature("docstring", cpt_loc_string) compute_local_real_energy;

%feature("autodoc", "fc(I, J, r)") fc;
%define fc_string
"
    Tapering function on pairwise contributions to total energy.

    Computes the species dependent tapering function on contributions to the total
    energy that depend on distances between pairs of atoms. Used in eq (3) - (6)
    of Zipoli & Curioni, New Journal of Physics 15 (2013) 123006.

    Parameters
    ----------

    I          : Species index of first atom (GE=0, TE=1).

    J          : Species index of second atom (GE=0, TE=1).

    r          : Distance between the two atoms (Angstrom).


"
%enddef
%feature("docstring", fc_string) fc;

%feature("autodoc", "fs(I, z)") fs;
%define fs_string
"
    Deviation from typical coordination.

    Computes the species-dependent deviation from typical coordination. 
    Defined in equation (10) of  Zipoli & Curioni, New Journal of 
    Physics 15 (2013) 123006.

    Parameters
    ----------

    I          : Species index of atom (GE=0, TE=1).

    z          : Coordination of atom.


"
%enddef
%feature("docstring", fs_string) fs;

%feature("autodoc", "tijk(ctheta, I, J, K)") tijk;
%define tijk_string
"
    Angular term.

    Computes the angular term associated with angle theta_{ijk} as defined
    in equation (11) of   Zipoli & Curioni, New Journal of Physics 15 (2013) 
    123006.

    Parameters
    ----------

    ctheta     : Cosine of angle theta_{ijk}

    I          : Species index of atom i (GE=0, TE=1).

    J          : Species index of atom j (GE=0, TE=1).

    K          : Species index of atom k (GE=0, TE=1).

"
%enddef
%feature("docstring", tijk_string) tijk;

%feature("autodoc", "tijk_fast(ctheta, I, J, K)") tijk_fast;
%define tijk_fast_string
"
    Angular term.

    Computes the angular term associated with angle theta_{ijk} as defined
    in equation (11) of   Zipoli & Curioni, New Journal of Physics 15 (2013) 
    123006. EXPERIMENTAL OPTIMISED VERSION FOR FAST COMPUTATION.

    Parameters
    ----------

    ctheta     : Cosine of angle theta_{ijk}

    I          : Species index of atom i (GE=0, TE=1).

    J          : Species index of atom j (GE=0, TE=1).

    K          : Species index of atom k (GE=0, TE=1).

"
%enddef
%feature("docstring", tijk_fast_string) tijk_fast;


/* This will be parsed to generate the wrapper */
%include "energy.h"

