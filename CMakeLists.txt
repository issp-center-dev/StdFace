cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(StdFace C)

option(UHF "Build uhf_dry.out" OFF)
MESSAGE(STATUS "Build uhf_dry.out - ${UHF}")

option(MVMC "Build mvmc_dry.out" OFF)
MESSAGE(STATUS "Build mvmc_dry.out - ${MVMC}")

option(HPHI "Build hphi_dry.out" OFF)
MESSAGE(STATUS "Build hphi_dry.out - ${HPHI}")

option(HWAVE "Build hwave_dry.out" OFF)
MESSAGE(STATUS "Build hwave_dry.out - ${HWAVE}")

add_subdirectory(src)
install(DIRECTORY samples DESTINATION share/stdface)


option(TestStdFace "Enable testing" ON)
if (TestStdFace)
  enable_testing()
  add_subdirectory(test)
endif()
