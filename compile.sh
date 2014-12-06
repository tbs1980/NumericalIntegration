mkdir -p test
cd test
mkdir -p testOutput
cd ..
mkdir -p build
cd build
cmake ..
cd ..
make -j7 -C build
