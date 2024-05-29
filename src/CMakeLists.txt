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
   add_executable(uhf_dry.out dry.c ${SOURCES_StdFace})
   target_compile_definitions(uhf_dry.out PUBLIC _UHF)
   target_link_libraries(uhf_dry.out PUBLIC m)
   install(TARGETS uhf_dry.out RUNTIME DESTINATION bin)
endif()

if(MVMC)
   add_executable(mvmc_dry.out dry.c ${SOURCES_StdFace})
   target_compile_definitions(mvmc_dry.out PUBLIC _mVMC)
   target_link_libraries(mvmc_dry.out PUBLIC m)
   install(TARGETS mvmc_dry.out RUNTIME DESTINATION bin)
endif()

if(HPHI)
   add_executable(hphi_dry.out dry.c ${SOURCES_StdFace})
   target_compile_definitions(hphi_dry.out PUBLIC _HPhi)
   target_link_libraries(hphi_dry.out PUBLIC m)
   install(TARGETS hphi_dry.out RUNTIME DESTINATION bin)
endif()

if (HWAVE)
   add_executable(hwave_dry.out dry.c ${SOURCES_StdFace})
   target_compile_definitions(hwave_dry.out PUBLIC _HWAVE)
   target_link_libraries(hwave_dry.out PUBLIC m)
   install(TARGETS hwave_dry.out RUNTIME DESTINATION bin)
endif()

add_library(StdFace STATIC ${SOURCES_StdFace})
