# - Find LAPACK library
# This module finds an installed fortran library that implements the LAPACK
# linear-algebra interface (see http://www.netlib.org/lapack/).
#
# The approach follows that taken for the autoconf macro file, acx_lapack.m4
# (distributed at http://ac-archive.sourceforge.net/ac-archive/acx_lapack.html).
#
# This module sets the following variables:
#  LAPACK_FOUND - set to true if a library implementing the LAPACK interface
#    is found
#  LAPACK_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#    and -L).
#  LAPACK_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use LAPACK
#  LAPACK95_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use LAPACK95
#  LAPACK95_FOUND - set to true if a library implementing the LAPACK f95
#    interface is found
#  BLAS_STATIC  if set on this determines what kind of linkage we do (static)
#  BLAS_VENDOR  if set checks only the specified vendor, if not set checks
#     all the possibilities
#  BLA_F95     if set on tries to find the f95 interfaces for BLAS/LAPACK
### List of vendors (BLA_VENDOR) valid in this module
##  Intel(mkl), ACML,Apple, NAS, Generic

#=============================================================================
# Copyright 2007-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

set(_lapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT _LANGUAGES_ MATCHES Fortran)
include( CheckFunctionExists )
else ()
include( CheckFortranFunctionExists)
endif ()

set(LAPACK_FOUND FALSE)
set(LAPACK95_FOUND FALSE)

# TODO: move this stuff to separate module

macro(Check_Lapack_Libraries LIBRARIES _prefix _name _flags _list _blas _threads)
  # This macro checks for the existence of the combination of fortran libraries
  # given by _list.  If the combination is found, this macro checks (using the
  # Check_Fortran_Function_Exists macro) whether can link against that library
  # combination using the name of a routine given by _name using the linker
  # flags given by _flags.  If the combination of libraries is found and passes
  # the link test, LIBRARIES is set to the list of complete library paths that
  # have been found.  Otherwise, LIBRARIES is set to FALSE.
  
  # N.B. _prefix is the prefix applied to the names of all cached variables that
  # are generated internally and marked advanced by this macro.
  
  set(_libraries_work TRUE)
  set(${LIBRARIES})
  set(_combined_name)

  if (NOT _libdir)
    if (WIN32)
      set(_libdir $ENV{LIB})
    elseif (APPLE)
      set(_libdir $ENV{DYLD_LIBRARY_PATH} /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64)
    else ()
      set(_libdir $ENV{LD_LIBRARY_PATH} /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64)
    endif ()
  endif ()

  foreach(_library ${_list})
    set(_combined_name ${_combined_name}_${_library})
    
    if(_libraries_work)
      if (BLAS_STATIC)
	if (WIN32)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
	endif ()
	if (APPLE)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
	else ()
          set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
	endif ()
      else ()
	if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
          # for ubuntu's libblas3gf and liblapack3gf packages
          set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
	endif ()
      endif ()
      find_library(${_prefix}_${_library}_LIBRARY
       	NAMES ${_library}
       	PATHS ${_libdir}
       	)
      #find_library(${_prefix}_${_library}_LIBRARY ${_library})
      mark_as_advanced(${_prefix}_${_library}_LIBRARY)
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
    endif()
  endforeach()


  
  if(_libraries_work)
    # Test this combination of libraries.
    if(UNIX AND BLAS_STATIC)
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} "-Wl,--start-group" ${${LIBRARIES}} ${_blas} "-Wl,--end-group" ${_threads})
    else()
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blas} ${_threads})
    endif()
    #  message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
    if (NOT _LANGUAGES_ MATCHES Fortran)
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    else ()
      check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
    endif ()
    set(CMAKE_REQUIRED_LIBRARIES)
    mark_as_advanced(${_prefix}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
    #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
  endif()
  
  if(_libraries_work)
    set(${LIBRARIES} ${${LIBRARIES}} ${_blas} ${_threads})
  else()
    set(${LIBRARIES} FALSE)
  endif()
  
endmacro()


set(LAPACK_LINKER_FLAGS)
set(LAPACK_LIBRARIES)
set(LAPACK95_LIBRARIES)


if(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
  find_package(BLAS)
else()
  find_package(BLAS REQUIRED)
endif()


if(BLAS_FOUND)
  set(LAPACK_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
  
  # GOTO
  if (BLAS_VENDOR_FOUND STREQUAL "GOTO")
    if(NOT LAPACK_LIBRARIES)
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACK
	cheev
	""
	"goto2"
	"${BLAS_LIBRARIES}"
	""
	)
    endif()
  endif ()

  # ACML
  if (BLAS_VENDOR_FOUND STREQUAL "ACML" OR BLAS_VENDOR_FOUND STREQUAL "ACML_MP")
    set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
  endif ()

  # ACCELERATE (APPLE)
  if (BLAS_VENDOR_FOUND STREQUAL "ACCELERATE")
    if(NOT LAPACK_LIBRARIES)
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACK
	cheev
	""
	"Accelerate"
	"${BLAS_LIBRARIES}"
	""
	)
    endif()
  endif ()

  # MKL
  if (BLAS_VENDOR_FOUND MATCHES "MKL")
    set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
  endif()

  # ATLAS
  if (BLAS_VENDOR_FOUND MATCHES "ATLAS")
    check_lapack_libraries(
      LAPACK_LIBRARIES
      LAPACK
      cheev
      ""
      "lapack;lapack_atlas"
      "${BLAS_LIBRARIES}"
      ""
      )
  endif()

  # GENERIC LAPACK IMPLEMENTATION
  if ( NOT LAPACK_LIBRARIES )
    check_lapack_libraries(
      LAPACK_LIBRARIES
      LAPACK
      cheev
      ""
      "lapack"
      "${BLAS_LIBRARIES}"
      ""
      )
  endif ()

else()
  message(STATUS "LAPACK requires BLAS")
endif()

if(LAPACK_LIBRARIES)
  set(LAPACK_FOUND TRUE)
else()
  set(LAPACK_FOUND FALSE)
endif()

if(NOT LAPACK_FIND_QUIETLY)
  if(LAPACK_FOUND)
    message(STATUS "A library with LAPACK API found.")
  else()
    if(LAPACK_FIND_REQUIRED)
      message(FATAL_ERROR
	"A required library with LAPACK API not found. Please specify library location."
	)
    else()
      message(STATUS
	"A library with LAPACK API not found. Please specify library location."
	)
    endif()
  endif()
endif()


set(CMAKE_FIND_LIBRARY_SUFFIXES ${_lapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
