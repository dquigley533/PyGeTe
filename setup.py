###########################
# Setup script for PyGeTe #
###########################
#!/usr/bin/env python

def configuration():

    from numpy.distutils.misc_util import Configuration
    import numpy

    config = Configuration('PyGeTe', parent_name=None, top_path=None)

    energy_src = ['src/util.f90','src/energy.f90']
    energy_inc = ['include/energy.h']

    config.add_library('energy', sources=energy_src)

    config.add_extension('_energy',
        sources = ['src/energy.i'],
        libraries = ['energy'],
        include_dirs = ['include'],
        depends = ['src/energy.i'] + energy_src + energy_inc,
    )

    config.version = "0.1.1"

    return config

if __name__ == "__main__":

    from numpy.distutils.core import setup


    setup(configuration=configuration,
          author       = "David Quigley",
          author_email = "d.quigley@warwick.ac.uk",
          description  = "Monte Carlo code for many body potential simulations of GeTe",
          url          = "https://github.com/dquigley-warwick/PyGeTe")



