rm -r build
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DGMX_PYTHON_PACKAGE=ON -DCMAKE_INSTALL_PREFIX=/home/oliverfl/git/gromacs/bin/
make
make check
sudo make install
source /home/oliverfl/git/gromacs/bin/bin/GMXRC
