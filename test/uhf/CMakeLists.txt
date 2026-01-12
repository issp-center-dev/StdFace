add_test(
  NAME HubbardSquare
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardSquare ${CMAKE_SOURCE_DIR}/test/uhf
)
add_test(
  NAME HubbardTriangle
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardTriangle ${CMAKE_SOURCE_DIR}/test/uhf
)
