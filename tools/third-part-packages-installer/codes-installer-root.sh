install_path_header=/opt/deepin15.6

# install mpich3
package=mpich3.23
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
#rm -rf ${packages[i]}


# install hdf5-mpich
package=mpich3.23
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


# install SuperLU-4.3
package=mpich3.23
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