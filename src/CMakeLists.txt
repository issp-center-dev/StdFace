# include guard
cmake_minimum_required(VERSION 2.8.0 )

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of StdFace.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)
include_directories(common)
set(SOURCES_StdFace
        StdFacemain.c output.c cal_energy.c green.c makeham.c diag.c initial.c matrixlapack.c readdef.c common/setmemory.c
 )

include_directories(sfmt)
add_definitions(-DMEXP=19937)

add_library(StdFace STATIC ChainLattice.c HoneycombLattice.c SquareLattice.c StdFace_main.c StdFace_ModelUtil.c TriangularLattice.c Ladder.c Kagome.c Orthorhombic.c Pyrochlore.c Wannier90.c FCOrtho.c setmemory.c export_wannier90.c)

if (UHF)
   add_executable(uhf_dry.out dry.c)
   target_link_libraries(uhf_dry.out StdFace m)
   install(TARGETS uhf_dry.out RUNTIME DESTINATION bin)
endif(UHF)

if (MVMC)
   add_executable(mvmc_dry.out dry.c)
   target_link_libraries(mvmc_dry.out StdFace m)
   install(TARGETS mvmc_dry.out RUNTIME DESTINATION bin)
endif(MVMC)

if (HPHI)
   add_executable(hphi_dry.out dry.c)
   target_link_libraries(hphi_dry.out StdFace m)
   install(TARGETS hphi_dry.out RUNTIME DESTINATION bin)
endif(HPHI)

if (HWAVE)
   add_executable(hwave_dry.out dry.c)
   target_link_libraries(hwave_dry.out StdFace m)
   install(TARGETS hwave_dry.out RUNTIME DESTINATION bin)
endif(HWAVE)

