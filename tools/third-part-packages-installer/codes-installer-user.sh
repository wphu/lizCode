:' for centos 6, set the following environments in the file ~/.bashrc
export PATH=/home/wphu/opt-gcc/gcc/bin:$PATH
export LD_LIBRARY_PATH=${HOME}/opt-gcc/gcc/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/opt-gcc/gmp/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/opt-gcc/mpc/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/opt-gcc/mpfr/lib:$LD_LIBRARY_PATH

before installing netcdf, first install: sudo apt install curl libcurl4-openssl-dev
'


#:' for gcc compilers
export install_path_header=${HOME}/opt-gcc
export compiler_c=gcc
export compiler_cxx=g++
export compiler_fortran=gfortran
export compiler_mpicc=mpicc
export compiler_mpicxx=mpicxx
export source_codes_root_path=$(pwd)
export compile_cores_number=20
#'


:' for intel compilers
source ${HOME}/opt-intel/intel/bin/compilervars.sh intel64

export install_path_header=${HOME}/opt-intel
export compiler_c=icc
export compiler_cxx=icc
export compiler_fortran=ifort
export compiler_mpicc=mpicc
export compiler_mpicxx=mpicxx
export source_codes_root_path=$(pwd)
export compile_cores_number=10
'


# install anaconda3
if [ -d ${install_path_header}/anaconda3 ];then
    echo "anaconda3 has been installed"
else
    bash ./Anaconda3-2018.12-Linux-x86_64.sh -b -p ${install_path_header}/anaconda3
fi


# install mpich3
export CC=${compiler_c}
export CXX=${compiler_cxx}
export FC=${compiler_fortran}
package=mpich-3.2.1
install_path=mpich
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
fi
export PATH=${install_path_header}/${install_path}/bin:$PATH



# install hdf5
export CC=${compiler_c}
export CXX=${compiler_cxx}
export FC=${compiler_fortran}
export CFLAGS=-fPIC
package=hdf5-1.8.20
install_path=hdf5
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
fi
export CFLAGS=""

# install hdf5-mpich
export CC=${compiler_mpicc}
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



# install fftw
export CC=${compiler_c}
package=fftw-3.3.4
install_path=fftw
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
fi


# install netcdf(netcdf-c)
export CC=${compiler_c}
export CPPFLAGS="-I${install_path_header}/hdf5/include"
export LDFLAGS="-L${install_path_header}/hdf5/lib"
package=netcdf-4.6.0
install_path=netcdf
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
    #rm -rf ${package}
fi
export CPPFLAGS=""
export LDFLAGS=""

# install netcdf-cxx
export CC=${compiler_c}
export CXX=${compiler_cxx}
export CPPFLAGS="-I${install_path_header}/hdf5/include -I${install_path_header}/netcdf/include"
export LDFLAGS="-L${install_path_header}/hdf5/lib -L${install_path_header}/netcdf/lib"
package=netcdf-cxx4-4.3.0
install_path=netcdf
if [ -d ${install_path_header}/${install_path} ];then
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
    #rm -rf ${package}
fi
export CPPFLAGS=""
export LDFLAGS=""

# install netcdf-cxx
export CC=${compiler_c}
export CPPFLAGS="-I${install_path_header}/hdf5/include -I${install_path_header}/netcdf/include"
export LDFLAGS="-L${install_path_header}/hdf5/lib -L${install_path_header}/netcdf/lib"
package=netcdf-cxx4-4.3.0
install_path=netcdf
if [ -d ${install_path_header}/${install_path} ];then
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
    #rm -rf ${package}
fi
export CPPFLAGS=""
export LDFLAGS=""

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

if [ -d ${install_path_header}/${install_path/lib64} ];then
    mv ${install_path_header}/${install_path}/lib64 ${install_path_header}/${install_path}/lib
fi

export FFLAGS=""
export OPTS=""
export DRVOPTS=""
export NOOPT=""

# install OpenBLAS
export FC=${compiler_fortran}
export F77=${compiler_fortran}
export FFLAGS=-fPIC
package=OpenBLAS-0.2.20
install_path=OpenBLAS
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
    cmake ../ -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ../..
fi
export FFLAGS=""

# install SuperLU
export CC=${compiler_c}
package=superlu_5.2.1
install_path=superlu
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    mv SuperLU_5.2.1 ${package}
    cd ${package}
    mkdir build
    cd build
    mkdir ${install_path_header}/${install_path}
    cmake ../ -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    if [ -d ${install_path_header}/${install_path/lib64} ];then
        mv ${install_path_header}/${install_path}/lib64 ${install_path_header}/${install_path}/lib
    fi
    cp CBLAS/libblas.a ${install_path_header}/${install_path}/lib/libblas.a
    cd ../..
fi



# install SuperLU-DIST
export CC=${compiler_mpicc}
package=superlu_dist_5.3.0
install_path=superlu_dist
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    mv SuperLU_DIST_5.3.0 ${package}
    cd ${package}
    mkdir build
    cd build
    mkdir ${install_path_header}/${install_path}
    cmake ../ \
    -DCMAKE_C_FLAGS="-std=c99 -g" \
    -Denable_blaslib=OFF \
    -Denable_parmetislib=OFF \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_C_COMPILER=${compiler_mpicc} \
    -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ../..
fi


# install PETSc
export CC=${compiler_c}
export FC=${compiler_fortran}
package=petsc-3.8.3
install_path=petsc
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    python2 ./configure --with-mpi-dir=${install_path_header}/mpich  --download-fblaslapack --prefix=${install_path_header}/${install_path}
    make  MAKE_NP=${compile_cores_number} all test
    make  install
    cd ..
fi


# sundials
export CC=${compiler_c}
export FC=${compiler_fortran}
package=sundials-2.7.0
install_path=sundials
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    mkdir sundials-build
    cd sundials-build

    cmake \
    -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path} \
    -DEXAMPLES_INSTALL_PATH=${install_path_header}/${install_path}/examples \
    -DCMAKE_LINKER=${install_path_header}/${install_path}/lib \
    -DLAPACK_ENABLE=ON \
    -DLAPACK_LIBRARIES="${install_path_header}/lapack/lib/libblas.a;${install_path_header}/lapack/lib/liblapack.a" \
    -DOPENMP_ENABLE=ON \
    -DMPI_ENABLE=ON \
    ../

    make -j${compile_cores_number}
    make install
    cd ../..
fi


# install umfpack included in SuiteSparse
export CC=${compiler_c}
export FC=${compiler_fortran}
export CFLAGS="-fPIC"
export CPPFLAGS="-fPIC"
package=SuiteSparse-5.3.0
install_path=SuiteSparse
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${source_codes_root_path}/${package}/lib"
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    mv SuiteSparse ${package}
    rm -rf SuiteSparse
    cd ${package}
    make  -j${compile_cores_number} BLAS="${install_path_header}/lapack/lib/libblas.a -lgfortran" LAPACK=${install_path_header}/lapack/lib/liblapack.a
    #make BLAS="/home/huwanpeng/source-codes/lapack-3.8.0/librefblas.a -lgfortran" LAPACK=/home/huwanpeng/source-codes/lapack-3.8.0/liblapack.a
    echo $LD_LIBRARY_PATH
    mkdir ${install_path_header}/${install_path}
    cp -r bin ${install_path_header}/${install_path}/bin
    cp -r lib ${install_path_header}/${install_path}/lib
    cp -r include ${install_path_header}/${install_path}/include
    cd ..
fi
export CFLAGS=""
export CPPFLAGS=""

# install gperftools
export CC=${compiler_c}
export CXX=${compiler_cxx}
package=gperftools-2.7
install_path=gperftools
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
fi