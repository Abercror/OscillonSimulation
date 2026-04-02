project(
    'OscillonSimulation',
    'cpp',
    version: '0.1.0',
    default_options: ['cpp_std=c++20'],
)

#opencl_inc = include_directories('/opt/homebrew/opt/opencl-clhpp-headers/include')
openmp_inc = include_directories('/opt/homebrew/opt/libomp/include')
hdf5_inc = include_directories('/opt/homebrew/opt/hdf5/include')
openmp_dep = dependency('openmp', required: true)
homebrew_inc = include_directories('/opt/homebrew/include')
include = include_directories('include')

subdir('src')