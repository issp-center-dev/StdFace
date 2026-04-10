add_test(
  NAME HubbardSquare
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardSquare ${CMAKE_SOURCE_DIR}/test/uhf
)
add_test(
  NAME HubbardTriangle
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardTriangle ${CMAKE_SOURCE_DIR}/test/uhf
)

# test python version
add_test(
  NAME python_HubbardSquare
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardSquare ${CMAKE_SOURCE_DIR}/test/uhf python
)
add_test(
  NAME python_HubbardTriangle
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardTriangle ${CMAKE_SOURCE_DIR}/test/uhf python
)
