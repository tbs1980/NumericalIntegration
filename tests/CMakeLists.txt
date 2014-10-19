cmake_minimum_required(VERSION 2.8)

project(NumericalIntegration_weights CXX)

if(UNIX)
	if(CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -Wfatal-errors -g")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
	endif(CMAKE_COMPILER_IS_GNUCXX)
	#is the compiler clang?
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -g")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++11 -pedantic -w3 -g")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=gnu++11 -O3 -march=native ")
	endif()
#############################
# NEED CHECKS FOR WINDOWS
#############################
endif(UNIX)

#Eigen3 headers
find_package(Eigen3 REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MPFRCPP REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${GMP_INCLUDES})
include_directories(${MPFR_INCLUDES})
include_directories(${MPFRCPP_INCLUDES})

# our headers
include_directories(${CMAKE_SOURCE_DIR})

# libs
set(ALL_LIBS ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

add_executable(testLaurieGautschi testLaurieGautschi.cpp)
add_executable(testQuadratureKronrod testQuadratureKronrod.cpp)
add_executable(testIntegrator_LogPow testIntegrator_LogPow.cpp)
add_executable(testIntegrator_Sine testIntegrator_Sine.cpp)
add_executable(testIntegrator_Peak testIntegrator_Peak.cpp)

target_link_libraries(testLaurieGautschi ${ALL_LIBS})
target_link_libraries(testQuadratureKronrod ${ALL_LIBS})
target_link_libraries(testIntegrator_LogPow ${ALL_LIBS})
target_link_libraries(testIntegrator_Sine ${ALL_LIBS})
target_link_libraries(testIntegrator_Peak ${ALL_LIBS})
