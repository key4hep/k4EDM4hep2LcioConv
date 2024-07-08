cmake_minimum_required(VERSION 3.14)
project(k4EDM4hep2LcioConv)

# project version
set( ${PROJECT_NAME}_VERSION_MAJOR 0 )
set( ${PROJECT_NAME}_VERSION_MINOR 8 )
set( ${PROJECT_NAME}_VERSION_PATCH 2 )
set( ${PROJECT_NAME}_VERSION  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}" )

# Define a default build type can be overriden by passing
# ``-DCMAKE_BUILD_TYPE=<type>`` when invoking CMake
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo
      CACHE STRING "Choose the type of build, options are: None Release MinSizeRel Debug RelWithDebInfo"
      FORCE
      )
  else()
    set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
      CACHE STRING "Choose the type of build, options are: None Release MinSizeRel Debug RelWithDebInfo"
      FORCE
      )
  endif()
endif()

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
include(cmake/build_flags.cmake)
k4edm4hep2lcioconv_set_compiler_flags()

find_package(LCIO 2.20.1 REQUIRED)
find_package(podio REQUIRED)
find_package(EDM4HEP 0.10.1 REQUIRED)
find_package(ROOT REQUIRED COMPONENTS MathCore)

add_subdirectory(k4EDM4hep2LcioConv)
add_subdirectory(standalone)

include(CTest)
if (BUILD_TESTING)
  add_subdirectory(tests)
endif()

include(cmake/k4EDM4hep2LcioConvCreateConfig.cmake)
