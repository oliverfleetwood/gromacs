rm -r build
mkdir build
cd build
export GMX_BUILD_UNITTESTS=OFF
export BUILD_TESTING=OFF
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DBUILD_TESTING=ON -DGMX_BUILD_UNITTESTS=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DGMX_PYTHON_PACKAGE=ON -DCMAKE_INSTALL_PREFIX=/home/oliverfl/git/gromacs/bin/ #-DCMAKE_BUILD_TYPE=debug
make
make check
sudo make install
source /home/oliverfl/git/gromacs/bin/bin/GMXRC


