# Python package container wrapped around PyGe
# library implementing the Zipoli many body
# potential.

author = 'David Quigley'

# Define some constants to refer to atom types
# as they are defined inside the Fortran lib
GE = 0
TE = 1

symdict = {"Ge":GE, "Te":TE}

from .energy import *
