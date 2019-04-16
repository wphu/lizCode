:' for centos 6, set the following environments in the file ~/.bashrc
export PATH=/home/wphu/opt-gcc/gcc/bin:$PATH
export LD_LIBRARY_PATH=/home/wphu/opt-gcc/gcc/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/wphu/opt-gcc/gmp/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/wphu/opt-gcc/mpc/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/wphu/opt-gcc/mpfr/lib:$LD_LIBRARY_PATH
'


#:' for gcc compilers
export install_path_header=${HOME}/opt-gcc
export compiler_c=gcc
export compiler_cxx=g++
export compiler_fortran=gfortran
export compiler_mpicc=mpicc
export compiler_mpicxx=mpicxx
export source_codes_root_path=$(pwd)
export compile_cores_number=10
#'


:' for intel compilers
source /home/huwanpeng/opt-intel/intel/bin/compilervars.sh intel64

export install_path_header=/home/huwanpeng/opt-intel
export compiler_c=icc
export compiler_cxx=icc
export compiler_fortran=ifort
export compiler_mpicc=mpicc
export compiler_mpicxx=mpicxx
export source_codes_root_path=$(pwd)
export compile_cores_number=10
'


# install lapack
export FORTRAN=${compiler_fortran}
export FFLAGS="-fPIC"
export OPTS="-O2 -frecursive"
export DRVOPTS=${OPTS}
export NOOPT="-O0 -frecursive"
package=lapack-3.8.0
install_path=lapack
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    mkdir build
    cd build
    mkdir ${install_path_header}/${install_path}
    cmake ../ -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path} -DBUILD_SHARED_LIBS=ON
    make -j${compile_cores_number}
    make install
    cd ../..
fi

if [ -d ${install_path_header}/${install_path}/lib64} ];then
    mv ${install_path_header}/${install_path}/lib64 ${install_path_header}/${install_path}/lib
fi

export FFLAGS=""
export OPTS=""
export DRVOPTS=""
export NOOPT=""