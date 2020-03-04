###########################
# Setup script for PyGeTe #
###########################
#!/usr/bin/env python

def configuration():

    from numpy.distutils.misc_util import Configuration
    import numpy

    config = Configuration('PyGeTe', parent_name=None, top_path=None)

    energy_src = ['src/util.f90','src/energy.f90']
    config.add_library('energy', sources=energy_src)

    config.add_extension('_energy',
        sources = ['src/energy.i'],
        libraries = ['energy'],
        include_dirs = ['include'],
        depends = ['src/energy.i']
    )

    config.version = "0.1.1"

    return config

if __name__ == "__main__":

    from numpy.distutils.core import setup

    setup(configuration=configuration, zip_safe=True)



