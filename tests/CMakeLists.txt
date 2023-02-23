add_library(edmCompare SHARED CompareEDM4hepLCIO.cc)
target_link_libraries(edmCompare PUBLIC EDM4HEP::edm4hep ${LCIO_LIBRARIES})
target_include_directories(edmCompare PUBLIC ${LCIO_INCLUDE_DIRS})
