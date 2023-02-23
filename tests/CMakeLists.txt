add_library(edmCompare SHARED src/CompareEDM4hepLCIO.cc)
target_link_libraries(edmCompare PUBLIC EDM4HEP::edm4hep ${LCIO_LIBRARIES})
target_include_directories(edmCompare PUBLIC ${LCIO_INCLUDE_DIRS})

add_executable(compare-contents compare_contents.cpp)
target_link_libraries(compare-contents PRIVATE edmCompare podio::podioRootIO)
target_include_directories(compare-contents PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>)
