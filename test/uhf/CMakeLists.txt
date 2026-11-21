add_test(
  NAME HubbardSquare
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardSquare ${CMAKE_SOURCE_DIR}/test/uhf
)
add_test(
  NAME HubbardTriangle
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardTriangle ${CMAKE_SOURCE_DIR}/test/uhf
)
add_test(
  NAME SpinChain
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh SpinChain ${CMAKE_SOURCE_DIR}/test/uhf
)
add_test(
  NAME SpinGCChain
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh SpinGCChain ${CMAKE_SOURCE_DIR}/test/uhf
)
