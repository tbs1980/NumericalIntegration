find_path(MPFRCPP_INCLUDES
  NAMES
  mpreal.h
  PATHS
  ${MPFRCPP_ROOT}
  ${INCLUDE_INSTALL_DIR}
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFRCPP DEFAULT_MSG MPFRCPP_INCLUDES)
mark_as_advanced(MPFRCPP_INCLUDES)
