add_library(k4EDM4hep2LcioConv SHARED
  src/k4EDM4hep2LcioConv.cpp)

target_include_directories(k4EDM4hep2LcioConv PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(k4EDM4hep2LcioConv PUBLIC
  ${LCIO_LIBRARIES}
  EDM4HEP::edm4hep)
