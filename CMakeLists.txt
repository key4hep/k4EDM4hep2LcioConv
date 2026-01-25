cmake_minimum_required(VERSION 3.14)
project(k4EDM4hep2LcioConv LANGUAGES CXX)

# project version
set( ${PROJECT_NAME}_VERSION_MAJOR 0 )
set( ${PROJECT_NAME}_VERSION_MINOR 14 )
set( ${PROJECT_NAME}_VERSION_PATCH 0 )
set( ${PROJECT_NAME}_VERSION  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}" )

include(cmake/Key4hepConfig.cmake)

include(GNUInstallDirs)

include(cmake/build_flags.cmake)
k4edm4hep2lcioconv_set_linker_flags()
find_package(LCIO 2.22 REQUIRED)
find_package(podio 1.0 REQUIRED)
find_package(EDM4HEP 1.0)
if (NOT EDM4HEP_FOUND)
  find_package(EDM4HEP 0.99 REQUIRED)
endif()
find_package(ROOT REQUIRED COMPONENTS MathCore)

add_subdirectory(k4EDM4hep2LcioConv)
add_subdirectory(standalone)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()

include(cmake/k4EDM4hep2LcioConvCreateConfig.cmake)
