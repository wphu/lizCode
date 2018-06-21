install_path_header=/opt/deepin15.6
compiler_c=gcc
compiler_fortran=gfortran

# install mpich3
CC=compiler_c
package=mpich-3.2.1
install_path=mpich
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
fi
tar -xvf ${package}.tar.gz
tar -xvf ${package}.tar
tar -xvf ${package}.gz
cd ${package}
./configure --prefix=${install_path_header}/${install_path}
make
make install
cd ..
export PATH=${install_path_header}/${install_path}:$PATH
which mpicc
#rm -rf ${packages[i]}


# install hdf5-mpich
CC=mpicc
package=hdf5-1.8.20
install_path=hdf5-mpich
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
fi
tar -xvf ${package}.tar.gz
tar -xvf ${package}.tar
tar -xvf ${package}.gz
cd ${package}
./configure --enable-parallel --prefix=${install_path_header}/${install_path}
make
make install
cd ..


# install SuperLU-4.3
CC=compiler_c
FC=compiler_fortran
package=superlu_4.3
install_path=SuperLU_4.3
if [ -d ${install_path_header}/${install_path} ];then
    echo "${package} has been installed"
    continue
fi
tar -xvf ${package}.tar.gz
tar -xvf ${package}.tar
tar -xvf ${package}.gz
cd ${package}
./configure --prefix=${install_path_header}/${install_path}
make
make install
cd ..
chmod a+r -R ${install_path_header}/${install_path