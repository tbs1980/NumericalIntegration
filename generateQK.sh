cd build
make generateGaussKronrodQuadrature
cd ..
./build/bin/generateGaussKronrodQuadrature
cp ./test/testOutput/GaussKronrodNodesWeights.h ./GaussKronrodNodesWeights.h

