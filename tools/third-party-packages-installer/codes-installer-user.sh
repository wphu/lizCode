packages=(mpich-3.2.1 hdf5-1.8.20)
username=wphu
install_path_header=/opt
install_paths=(mpich-intel hdf5-intel)



packages_number=${#packages[*]}
for i in $(seq 0 `expr $packages_number - 1`)
do
    if [ -d /home/${username}${install_path_header}/${install_paths[i]} ];then
        echo "${packages[i]} has been installed"
        continue
    fi

    tar -xvf ${packages[i]}.tar.gz
    tar -xvf ${packages[i]}.tar
    tar -xvf ${packages[i]}.gz
    cd ${packages[i]}
    ./configure --prefix=/home/${username}${install_path_header}/${install_paths[i]}
    make
    make install
    cd ..
    rm -rf ${packages[i]}
done

