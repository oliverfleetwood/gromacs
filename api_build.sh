source activate py37
python -m ensurepip --default-pip
pip install --upgrade pip setuptools
pip install --upgrade cmake scikit-build
source /home/oliverfl/git/gromacs/bin/bin/GMXRC
pip install gmxapi