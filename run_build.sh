rm -rf build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$(realpath ../install) ..
make -j 5 install
