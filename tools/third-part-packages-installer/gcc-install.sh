// before executing the script, execute the following command
// sudo apt install gcc g++ gfortran m4 make cmake


export install_path_header=/home/huwanpeng/opt-gcc
export compiler_c=gcc
export compiler_cxx=g++
export compiler_fortran=gfortran
export compiler_mpicc=mpicc
export compiler_mpicxx=mpicxx
export source_codes_root_path=$(pwd)
export compile_cores_number=10

# install gmp
export CC=${compiler_c}
#export CPPFLAGS="-I${install_path_header}/hdf5/include -I${install_path_header}/netcdf/include"
#export LDFLAGS="-L${install_path_header}/hdf5/lib -L${install_path_header}/netcdf/lib"
package=gmp-6.1.2
install_path=gmp
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    tar -xvf ${package}.tar.lz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path}
    make -j${compile_cores_number}
    make install
    cd ..
fi


# install mpfr
export CC=${compiler_c}
#export CPPFLAGS="-I${install_path_header}/hdf5/include -I${install_path_header}/netcdf/include"
#export LDFLAGS="-L${install_path_header}/hdf5/lib -L${install_path_header}/netcdf/lib"
package=mpfr-4.0.1
install_path=mpfr
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    tar -xvf ${package}.tar.lz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path} --with-gmp=${install_path_header}/gmp
    make -j${compile_cores_number}
    make install
    cd ..
fi


# install mpc
export CC=${compiler_c}
#export CPPFLAGS="-I${install_path_header}/hdf5/include -I${install_path_header}/netcdf/include"
#export LDFLAGS="-L${install_path_header}/hdf5/lib -L${install_path_header}/netcdf/lib"
package=mpc-1.1.0
install_path=mpc
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    tar -xvf ${package}.tar.lz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path} --with-gmp=${install_path_header}/gmp --with-mpfr=${install_path_header}/mpfr
    make -j${compile_cores_number}
    make install
    cd ..
fi


# install gcc
export CC=${compiler_c}
#export CPPFLAGS="-I${install_path_header}/hdf5/include -I${install_path_header}/netcdf/include"
#export LDFLAGS="-L${install_path_header}/hdf5/lib -L${install_path_header}/netcdf/lib"
package=gcc-8.1.0
install_path=gcc
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
else
    tar -xvf ${package}.tar.gz
    tar -xvf ${package}.tar
    tar -xvf ${package}.gz
    tar -xvf ${package}.tar.lz
    cd ${package}
    ./configure --prefix=${install_path_header}/${install_path} --with-gmp=${install_path_header}/gmp --with-mpfr=${install_path_header}/mpfr --with-mpc=${install_path_header}/mpc --enable-threads=posix --disable-checking --enable--long-long --enable-languages=c,c++,fortran --disable-multilib
    make -j${compile_cores_number}
    make install
    cd ..
fi


export LD_LIBRARY_PATH=/home/huwanpeng/opt-gcc/gcc/lib/../lib64:$LD_LIBRARY_PATH
export LD_RUN_PATH=/home/huwanpeng/opt-gcc/gcc/lib/../lib64:$LD_RUN_PATH
export PATH=/home/huwanpeng/opt-gcc/gcc/bin:$PATH