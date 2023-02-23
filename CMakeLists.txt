cmake_minimum_required(VERSION 3.14)
project(k4EDM4hep2LcioConv)

# project version
set( ${PROJECT_NAME}_VERSION_MAJOR 0 )
set( ${PROJECT_NAME}_VERSION_MINOR 3 )
set( ${PROJECT_NAME}_VERSION_PATCH 0 )
set( ${PROJECT_NAME}_VERSION  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}" )

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

if(NOT CMAKE_CXX_STANDARD MATCHES "17|20")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

# Ninja compiler output
include(cmake/compiler_output.cmake)

find_package(LCIO REQUIRED)
find_package(EDM4HEP REQUIRED)

add_subdirectory(k4EDM4hep2LcioConv)
add_subdirectory(standalone)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()

include(cmake/k4EDM4hep2LcioConvCreateConfig.cmake)
