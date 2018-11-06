export install_path_header=/home/huwanpeng/opt
export compiler_c=gcc
export compiler_cxx=g++
export compiler_fortran=gfortran
export compiler_mpicc=mpicc
export compiler_mpicxx=mpicxx
export source_codes_root_path=$(pwd)

# install anaconda3
bash ./Anaconda3-5.1.0-Linux-x86_64.sh -b -p ${install_path_header}/anaconda3

# install mpich3
CC=${compiler_c}
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
    make
    make install
    cd ..
fi
export PATH=${install_path_header}/${install_path}/bin:$PATH



# install hdf5
CC=${compiler_c}
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
    make
    make install
    cd ..
fi
export CFLAGS=""

# install hdf5-mpich
CC=${compiler_mpicc}
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
    make
    make install
    cd ..
fi



# install fftw
CC=${compiler_c}
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
    make
    make install
    cd ..
fi


# install netcdf(netcdf-c)
CC=${compiler_c}
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
    make
    make install
    cd ..
    #rm -rf ${package}
fi
export CPPFLAGS=""
export LDFLAGS=""

# install netcdf-cxx
CC=${compiler_c}
CXX=${compiler_cxx}
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
    make
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
    make
    make install
    cd ..
    #rm -rf ${package}
fi
export CPPFLAGS=""
export LDFLAGS=""

# install lapack
FORTRAN=${compiler_fortran}
export FFLAGS="-fPIC"
OPTS="-O2 -frecursive"
DRVOPTS=${OPTS}
NOOPT="-O0 -frecursive"
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
    cmake ../ -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path}
    make
    make install
    cd ../..
fi
export FFLAGS=""
export OPTS=""
export DRVOPTS=""
export NOOPT=""

# install OpenBLAS
FC=${compiler_fortran}
F77=${compiler_fortran}
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
    make
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
    make
    make install
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
    make
    make install
    cd ../..
fi


# install PETSc
CC=${compiler_c}
FC=${compiler_fortran}
package=petsc-3.8.3
install_path=petsc
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    cd ${package}
    python2 ./configure --with-mpi-dir=${install_path_header}/mpich --download-fblaslapack --prefix=${install_path_header}/${install_path}
    make PETSC_DIR=${source_codes_root_path}/${package} PETSC_ARCH=arch-linux2-c-debug all
    make PETSC_DIR=${source_codes_root_path}/${package} PETSC_ARCH=arch-linux2-c-debug install
    cd ..
fi


# sundials
CC=${compiler_c}
FC=${compiler_fortran}
package=sundials-3.1.0
install_path=sundials
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    mkdir examples
    mkdir lib
    mkdir sundials-build

    cmake \
    -DCMAKE_INSTALL_PREFIX=${install_path_header}/${install_path} \
    -DEXAMPLES_INSTALL_PATH=${install_path_header}/${install_path}/examples \
    -DCMAKE_LINKER=${install_path_header}/${install_path}/lib \
    -DLAPACK_ENABLE=ON \
    -DOPENMP_ENABLE=ON \
    -DMPI_ENABLE=ON \
    ./${package}

    make
    make install
    rm -rf ${install_path_header}/examples ${install_path_header}/lib ${install_path_header}/sundials-build
fi


# install umfpack included in SuiteSparse
CC=${compiler_c}
FC=${compiler_fortran}
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
    make BLAS="${install_path_header}/lapack/lib/libblas.a -lgfortran" LAPACK=${install_path_header}/lapack/lib/liblapack.a
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

