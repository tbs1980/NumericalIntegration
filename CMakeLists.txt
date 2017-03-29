CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(NumericalIntegration CXX)

# Guard against in-source builds (got this from Eigen)
IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
ENDIF()

#################################################
# Compiler checks for Unix type systems.
#################################################
#IF(UNIX)
#	IF(CMAKE_COMPILER_IS_GNUCXX)
#		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -Wfatal-errors -g")
#		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
#		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
#	ENDIF(CMAKE_COMPILER_IS_GNUCXX)
#
#	#is the compiler clang?
#	IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -g")
#		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
#		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
#	ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++11 -pedantic -w3 -g")
#		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
#		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=gnu++11 -O3 -march=native ")
#	ENDIF()
#
#ENDIF(UNIX)

#################################################
# NEED CHECKS FOR WINDOWS
#################################################


INCLUDE(ExternalProject)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#Eigen3 headers
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(GMP REQUIRED)
FIND_PACKAGE(MPFR REQUIRED)
FIND_PACKAGE(MPFRCPP REQUIRED)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GMP_INCLUDES})
INCLUDE_DIRECTORIES(${MPFR_INCLUDES})
INCLUDE_DIRECTORIES(${MPFRCPP_INCLUDES})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# libraries
SET(ALL_LIBS ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

ADD_DEFINITIONS(
  -O3
  -DAPL=0
  -DIBM=0
  -DLIN=1
  -DNDEBUG
  -DCOMPUTE_NI_WEIGHTS_ON_THE_FLY

  -Wall
  -Wcast-align
  -Wchar-subscripts
  -Wconversion
  -Wdouble-promotion
  -Wextra
  -Wfatal-errors
  -Wformat-security
  -Wlogical-op
  -Wno-long-long  
  -Wno-psabi
  -Wno-variadic-macros
  #-Wnon-virtual-dtor
  -Wpointer-arith
  #-Wundef
  -Wuninitialized
  -Wunused-local-typedefs
  -Wwrite-strings

  -fstrict-aliasing
  -fno-check-new
  -fno-common
  -fvisibility=hidden
  -fPIC

  -c
  -g
  -march=native
  -pedantic
  -pthread
  )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LD_VERSION_SCRIPT_FLAG}")

ENABLE_TESTING()
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(test)
