#!/usr/bin/env bash
if [ ${TRAVIS_PYTHON_VERSION} = "2.7" ]; then
    sudo apt-get install -qq swig libeigen3-dev python-dev git-core python-numpy
elif [ ${TRAVIS_PYTHON_VERSION:0:1} = 3 ]; then
    sudo apt-get install -qq swig libeigen3-dev python-dev git-core python3-numpy python3-nose
else
    echo "ERROR: Unknown Python version. Got '${TRAVIS_PYTHON_VERSION}' from environment."
fi

git clone https://github.com/openbabel/openbabel.git openbabel_source

cd openbabel_source
git checkout ${HASH}
cd ..

mkdir openbabel_build
cd openbabel_build
cmake -DCMAKE_CXX_FLAGS="-march=native -mno-avx" -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON ../openbabel_source
make -j 2
sudo make install
export PYTHONPATH=/usr/local/lib/python${TRAVIS_PYTHON_VERSION}/site-packages:${PYTHONPATH}
