# NumericalIntegration

Comparison of Piessens' method to Laurie-Gautschi approach

## Requirements

	* [Eigen library](http://eigen.tuxfamily.org/index.php?title=Main_Page )
	* [The GNU Multiple Precision Arithmetic Library](https://gmplib.org/ )
	* [The GNU MPFR Library](http://www.mpfr.org/ )
	* [MPFR C++](http://www.holoborodko.com/pavel/mpfr/ )

## Compilation
A compile script has been added to the top level directory, simply run: ./compile.sh

The following compilation flags must be passed

	* -DEIGEN3_INCLUDE_DIR
	* -DGMP_ROOT
	* -DMPFR_ROOT
	* -DMPFRCPP_ROOT

	$ mkdir build
	$ cd build
	$ cmake -DEIGEN3_INCLUDE_DIR=path_to_Eigen3 -DGMP_ROOT=path_GMP_root_dir -DMPFR_ROOT=path_to_MPFR_root_dir -DMPFRCPP_ROOT=path_to_MPFRC++_root_dir  path_to_GaussKronrad
	$ make

For example,

	$ cmake ../ -DEIGEN3_INCLUDE_DIR=/arxiv/libraries/ubuntu/gcc/eigen-3.2.1/include/eigen3 -DGMP_ROOT=/arxiv/libraries/ubuntu/gcc/gmp-6.0.0 -DMPFR_ROOT=/arxiv/libraries/ubuntu/gcc/mpfr-3.1.2 -DMPFRCPP_ROOT=/arxiv/libraries/ubuntu/gcc/mpfrc++-3.5.9
	$ make
