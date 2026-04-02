src_files = [
    'main.cpp',
]

executable(
    'OscillonSimulation',
    src_files,
    dependencies: [openmp_dep],
    include_directories: [homebrew_inc, openmp_inc, hdf5_inc, include],
    install : true
)