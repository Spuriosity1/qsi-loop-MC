project('dilutedlat', ['cpp', 'c'],
  version: '1.0.0',
  meson_version : '>= 1.5.0',
  default_options : [ 'warning_level=3', 'buildtype=debug', 'cpp_std=c++2a' ]
)

if get_option('buildtype') == 'release'
  add_project_arguments('-Ofast',  language : 'cpp')
endif

extra_hdf5_libs = declare_dependency(
  link_args: ['-lz', '-ldl', '-lm']
)

latlib_dep = dependency('liblatindex', required: true)
fftw_dep = dependency('fftw3', version: '>= 3.0', required: true)
h5_dep = dependency('hdf5', required: true, method:'pkg-config')

loopMC = executable('loopMC', 
  files(['loopMC.cpp', 'MC_routines.cpp', 'observables.cpp']),
  dependencies: [extra_hdf5_libs, latlib_dep, fftw_dep, h5_dep],
  include_directories: 'include'
  )

