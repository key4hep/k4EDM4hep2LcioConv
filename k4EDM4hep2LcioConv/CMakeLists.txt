add_library(k4EDM4hep2LcioConv SHARED
  src/k4EDM4hep2LcioConv.cpp
  src/k4Lcio2EDM4hepConv.cpp
  )
add_library(k4EDM4hep2LcioConv::k4EDM4hep2LcioConv ALIAS k4EDM4hep2LcioConv)

target_include_directories(k4EDM4hep2LcioConv PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(k4EDM4hep2LcioConv PUBLIC
  LCIO::lcio
  EDM4HEP::edm4hep
  EDM4HEP::edm4hepUtils
)

set(public_headers
  include/${PROJECT_NAME}/k4EDM4hep2LcioConv.h
  include/${PROJECT_NAME}/k4EDM4hep2LcioConv.ipp
  include/${PROJECT_NAME}/k4Lcio2EDM4hepConv.h
  include/${PROJECT_NAME}/k4Lcio2EDM4hepConv.ipp
  include/${PROJECT_NAME}/MappingUtils.h
)

set_target_properties(${PROJECT_NAME}
  PROPERTIES
  PUBLIC_HEADER "${public_headers}"
  )

install(TARGETS k4EDM4hep2LcioConv
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
  COMPONENT dev)
