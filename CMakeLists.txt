cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(StdFace C)

option(URF "Build urf_dri.out" ON)
MESSAGE(STATUS "Build urf_dri.out - ${URF}")

option(MVMC "Build mvmc_dri.out" ON)
MESSAGE(STATUS "Build mvmc_dri.out - ${MVMC}")

option(HPHI "Build hphi_dri.out" ON)
MESSAGE(STATUS "Build hphi_dri.out - ${HPHI}")

option(HWAVE "Build hwave_dri.out" ON)
MESSAGE(STATUS "Build hwave_dri.out - ${HWAVE}")

add_subdirectory(src)
install(DIRECTORY samples DESTINATION share/stdface)
