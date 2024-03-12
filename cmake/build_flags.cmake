#---------------------------------------------------------------------------------------------------
#---k4edm4hep2lcioconv_set_compiler_flags
#
#  Set compiler and linker flags
#
#---------------------------------------------------------------------------------------------------

macro(k4edm4hep2lcioconv_set_compiler_flags)
  include(CheckCXXCompilerFlag)

  SET(COMPILER_FLAGS -fPIC -Wall -Wextra -Wpedantic)

  # AppleClang/Clang specific warning flags
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    set ( COMPILER_FLAGS ${COMPILER_FLAGS} -Winconsistent-missing-override -Wheader-hygiene )
  endif()

  FOREACH( FLAG ${COMPILER_FLAGS} )
    ## need to replace the minus or plus signs from the variables, because it is passed
    ## as a macro to g++ which causes a warning about no whitespace after macro
    ## definition
    STRING(REPLACE "-" "_" FLAG_WORD ${FLAG} )
    STRING(REPLACE "+" "P" FLAG_WORD ${FLAG_WORD} )

    CHECK_CXX_COMPILER_FLAG( "${FLAG}" CXX_FLAG_WORKS_${FLAG_WORD} )
    IF( ${CXX_FLAG_WORKS_${FLAG_WORD}} )
      message( STATUS "Adding ${FLAG} to CXX_FLAGS" )
      SET ( CMAKE_CXX_FLAGS "${FLAG} ${CMAKE_CXX_FLAGS} ")
    ELSE()
      message( STATUS "NOT Adding ${FLAG} to CXX_FLAGS" )
    ENDIF()
  ENDFOREACH()

  # resolve which linker we use
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -Wl,--version OUTPUT_VARIABLE stdout ERROR_QUIET)
  if("${stdout}" MATCHES "GNU ")
    set(LINKER_TYPE "GNU")
    message( STATUS "Detected GNU compatible linker" )
  else()
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -Wl,-v ERROR_VARIABLE stderr )
    if(("${stderr}" MATCHES "PROGRAM:ld") AND ("${stderr}" MATCHES "PROJECT:ld64"))
      set(LINKER_TYPE "APPLE")
      message( STATUS "Detected Apple linker" )
    else()
      set(LINKER_TYPE "unknown")
      message( STATUS "Detected unknown linker" )
    endif()
  endif()

  if("${LINKER_TYPE}" STREQUAL "APPLE")
    SET ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-undefined,dynamic_lookup")
  elseif("${LINKER_TYPE}" STREQUAL "GNU")
    SET ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--allow-shlib-undefined")
  else()
    MESSAGE( WARNING "No known linker (GNU or Apple) has been detected, pass no flags to linker" )
  endif()

endmacro(k4edm4hep2lcioconv_set_compiler_flags)
