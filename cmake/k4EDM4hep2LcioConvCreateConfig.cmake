# - Use CMake's module to help generating relocatable config files
include(CMakePackageConfigHelpers)

# - Versioning
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/k4EDM4hep2LcioConvConfigVersion.cmake
  VERSION k4EDM4hep2LcioConv_VERSION
  COMPATIBILITY SameMajorVersion)

# - Install time config and target files
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/k4EDM4hep2LcioConvConfig.cmake.in
  "${PROJECT_BINARY_DIR}/k4EDM4hep2LcioConvConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/k4EDM4hep2LcioConv"
  PATH_VARS
    CMAKE_INSTALL_BINDIR
    CMAKE_INSTALL_INCLUDEDIR
    CMAKE_INSTALL_LIBDIR
  )

# - install and export
install(FILES
  "${PROJECT_BINARY_DIR}/k4EDM4hep2LcioConvConfigVersion.cmake"
  "${PROJECT_BINARY_DIR}/k4EDM4hep2LcioConvConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/k4EDM4hep2LcioConv"
  )
install(EXPORT k4EDM4hep2LcioConvTargets
  NAMESPACE k4EDM4hep2LcioConv::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/k4EDM4hep2LcioConv"
  )

