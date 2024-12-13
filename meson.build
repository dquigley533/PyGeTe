# meson.build for PyGeTe
# D. Quigley - University of Warwick

project(
	'PyGeTe',         # The name of the project
	['fortran', 'c'], # The languages that Meson should initialize

	# Keyword arguments:
	version : '1.1.0',
	default_options : ['optimization=2', 'debug=false', 'warning_level=3'],
)

# List of options to pass to fortran compiler
# add_project_arguments(['-fbacktrace', '-fbounds-check'], language:'fortran') 

# Use swig to generate wrapper code
swig = find_program('swig')
incdir = include_directories('include')

wrapper_code = custom_target('wrapper_code',
	     input : ['src/energy.i'],
	     output : ['energy_wrap.c', 'energy.py'],
	     command : [swig, '-python', '-I../include/', '-outcurrentdir', '@INPUT@'],
	     )

# Need Python headers to compile the wrapper code
pydep = dependency('python3', version : '>=3.0.0')

py = import('python').find_installation(pure: false)

# Source files for fortran library
energy_src = ['src/util.f90','src/energy.f90']

# Build fortran source and wrapper code into extension module
py.extension_module(
  '_energy',
  [energy_src, wrapper_code],
  install: true,
  include_directories : incdir,
)

# Finally define and install python sources
python_sources = [
  '__init__.py',
]
py.install_sources(python_sources, subdir : 'PyGeTe')