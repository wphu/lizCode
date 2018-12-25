export install_path_header=/home/wphu/opt-intel2018
export compiler_c=icc
export compiler_cxx=icc
export compiler_fortran=ifort
export compiler_mpicc=mpicc
export compiler_mpicxx=mpicxx
export compiler_mpifortran=mpif90
export source_codes_root_path=$(pwd)
export compile_cores_number=10

# install hdf5-mpich
export CC=${compiler_mpicc}
export FC=${compiler_mpifortran}
package=hdf5-1.8.20
install_path=hdf5-mpich
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --enable-parallel --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
fi

