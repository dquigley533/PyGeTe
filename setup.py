import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def configuration():

    from numpy.distutils.misc_util import Configuration

    config = Configuration('PyGeTe', parent_name=None, top_path=None)

    zipoli_src = ['src/util.f90','src/energy.f90']
    config.add_library('zipoli', sources=zipoli_src)

    sources  = ['src/energy.i']
    includes = ['include/energy.h']

    config.add_extension('_energy',
        sources = sources,
        libraries = ['zipoli'],
        include_dirs = ['include'],
        depends = [sources, includes])

    return config

if __name__ == "__main__":

    from numpy.distutils.core import setup

    setup(configuration=configuration, zip_safe=True)



