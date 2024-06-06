# include guard
cmake_minimum_required(VERSION 3.0)

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of StdFace.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)
include_directories(common)
include_directories(sfmt)
add_definitions(-DMEXP=19937)

set(SOURCES_StdFace
  ChainLattice.c HoneycombLattice.c SquareLattice.c StdFace_main.c StdFace_ModelUtil.c TriangularLattice.c Ladder.c Kagome.c Orthorhombic.c Pyrochlore.c Wannier90.c FCOrtho.c setmemory.c export_wannier90.c
)

if(UHF)
   add_library(StdFace_uhf STATIC ${SOURCES_StdFace})
   target_compile_definitions(StdFace_uhf PUBLIC _UHF)
   add_executable(uhf_dry.out dry.c)
   target_compile_definitions(uhf_dry.out PUBLIC _UHF)
   target_link_libraries(uhf_dry.out PUBLIC StdFace_uhf m)
   install(TARGETS uhf_dry.out RUNTIME DESTINATION bin)
endif()

if(MVMC)
   add_library(StdFace_mvmc STATIC ${SOURCES_StdFace})
   target_compile_definitions(StdFace_mvmc PUBLIC _mVMC)
   add_executable(mvmc_dry.out dry.c)
   target_compile_definitions(mvmc_dry.out PUBLIC _mVMC)
   target_link_libraries(mvmc_dry.out PUBLIC StdFace_mvmc m)
   install(TARGETS mvmc_dry.out RUNTIME DESTINATION bin)
endif()

if(HPHI)
   add_library(StdFace_hphi STATIC ${SOURCES_StdFace})
   target_compile_definitions(StdFace_hphi PUBLIC _HPhi)
   add_executable(hphi_dry.out dry.c)
   target_compile_definitions(hphi_dry.out PUBLIC _HPhi)
   target_link_libraries(hphi_dry.out PUBLIC StdFace_hphi m)
   install(TARGETS hphi_dry.out RUNTIME DESTINATION bin)
endif()

if (HWAVE)
   add_library(StdFace_hwave STATIC ${SOURCES_StdFace})
   target_compile_definitions(StdFace_hwave PUBLIC _HWAVE)
   add_executable(hwave_dry.out dry.c)
   target_compile_definitions(hwave_dry.out PUBLIC _HWAVE)
   target_link_libraries(hwave_dry.out PUBLIC StdFace_hwave m)
   install(TARGETS hwave_dry.out RUNTIME DESTINATION bin)
endif()

