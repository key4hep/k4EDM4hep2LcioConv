add_executable(lcio2edm4hep lcio2edm4hep.cpp)
target_link_libraries(lcio2edm4hep PRIVATE k4EDM4hep2LcioConv podio::podioRootIO)

install(TARGETS lcio2edm4hep
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
