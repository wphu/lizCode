export install_path_header=/home/huwanpeng/opt-gcc
export compiler_c=gcc
export compiler_cxx=g++
export compiler_fortran=gfortran
export compiler_mpicc=mpicc
export compiler_mpicxx=mpicxx
export source_codes_root_path=$(pwd)


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
    make
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
    make
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
    make
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
    make
    make install
    cd ..
fi