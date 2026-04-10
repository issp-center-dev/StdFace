add_test(
  NAME uhfr_first_step
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh uhfr_first_step ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME uhfk_cdwsdw_input_1_cdw
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh uhfk_cdwsdw_input_1_cdw ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME uhfk_cdwsdw_input_1_sdw
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh uhfk_cdwsdw_input_1_sdw ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME uhfk_dos_plot
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh uhfk_dos_plot ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME uhfk_finite_t
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh uhfk_finite_t ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME uhfk_first_step
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh uhfk_first_step ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME rpa_kobayashi_2004_coulomb
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh rpa_kobayashi_2004_coulomb ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME rpa_sample
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh rpa_sample ${CMAKE_SOURCE_DIR}/test/hwave
)
