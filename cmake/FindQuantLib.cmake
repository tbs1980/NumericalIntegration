if (QuantLib_INCLUDES AND QuantLib_LIBRARIES)
  set(QuantLib_FIND_QUIETLY TRUE)
endif (QuantLib_INCLUDES AND QuantLib_LIBRARIES)

find_path(QuantLib_INCLUDES
  NAMES
  ql/version.hpp
  PATHS
  ${QuantLib_ROOT}/include
  ${INCLUDE_INSTALL_DIR}
)

find_library(QuantLib_LIBRARIES QuantLib PATHS $ENV{QuantLibDIR} ${QuantLib_ROOT}/lib ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QuantLib DEFAULT_MSG
                                  QuantLib_INCLUDES QuantLib_LIBRARIES)
mark_as_advanced(QuantLib_INCLUDES QuantLib_LIBRARIES)
