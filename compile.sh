cd test
mkdir -p testOutput
cd ..
mkdir -p build
cd build
cmake ..
cmake -G Sublime\ Text\ 2\ -\ Unix\ Makefiles .. 
cd ..
make -j7 -C build
