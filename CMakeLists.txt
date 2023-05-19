cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(StdFace NONE)

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
enable_language(C)

if(UHF)
add_definitions(-D_UHF)
endif(UHF)

if(MVMC)
add_definitions(-D_mVMC)
endif(MVMC)

if(HPHI)
add_definitions(-D_HPhi)
endif(HPHI)

if(HWAVE)
add_definitions(-D_HWAVE)
endif(HWAVE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
endif()

find_package(LAPACK)
if(USE_SCALAPACK MATCHES OFF)
  if(LAPACK_FOUND)
    add_definitions(-D_lapack)
  endif(LAPACK_FOUND)
endif()
add_subdirectory(src)

# Build and enable tests
# testing setup
# enable_testing() must be called in the top-level CMakeLists.txt before any add_subdirectory() is called.
option(Testing "Enable testing" OFF)
