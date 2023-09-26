add_library(edmCompare SHARED src/CompareEDM4hepLCIO.cc src/ObjectMapping.cc)
target_link_libraries(edmCompare PUBLIC EDM4HEP::edm4hep ${LCIO_LIBRARIES})
target_include_directories(edmCompare PUBLIC ${LCIO_INCLUDE_DIRS})

add_library(TestUtils SHARED src/EDM4hep2LCIOUtilities.cc)
target_link_libraries(TestUtils PUBLIC EDM4HEP::edm4hep ${LCIO_LIBRARIES})
target_include_directories(TestUtils PUBLIC ${LCIO_INCLUDE_DIRS})

add_executable(compare-contents compare_contents.cpp)
target_link_libraries(compare-contents PRIVATE edmCompare podio::podioRootIO)
target_include_directories(compare-contents PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>)

add_executable(edm4hep_to_lcio edm4hep_to_lcio.cpp)
target_link_libraries(edm4hep_to_lcio PRIVATE edmCompare TestUtils k4EDM4hep2LcioConv)
target_include_directories(edm4hep_to_lcio PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>)

add_test(NAME edm4hep_to_lcio COMMAND edm4hep_to_lcio)

find_program(BASH_PROGRAM bash)

add_test(fetch_test_inputs ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/get_test_data.sh)

add_test(standalone_ild_rec_file ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_standalone_converter.sh ild_higgs_rec.slcio)

add_test(standalone_ild_dst_file ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_standalone_converter.sh ild_higgs_dst.slcio)

set_tests_properties(
    fetch_test_inputs
    standalone_ild_rec_file
    standalone_ild_dst_file
  PROPERTIES
    ENVIRONMENT "TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR};PATH=${CMAKE_CURRENT_BINARY_DIR}:${PROJECT_BINARY_DIR}/standalone:$ENV{PATH}"
)

set_tests_properties(
    standalone_ild_rec_file
    standalone_ild_dst_file
  PROPERTIES
    DEPENDS fetch_test_inputs
)
