add_test(
  NAME HubbardSquare
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardSquare ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME HubbardTriangle
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh HubbardTriangle ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)

# --- matrix tests (auto-generated) ---
add_test(
  NAME matrix_chain_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_chain_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_chain_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_chain_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_chain_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_chain_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_chain_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_chain_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_chain_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_chain_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_chain_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_chain_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_square_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_square_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_square_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_square_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_square_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_square_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_square_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_square_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_square_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_square_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_square_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_square_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_triangular_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_triangular_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_triangular_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_triangular_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_triangular_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_triangular_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_triangular_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_triangular_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_triangular_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_triangular_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_triangular_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_triangular_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_honeycomb_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_honeycomb_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_honeycomb_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_honeycomb_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_honeycomb_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_honeycomb_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_honeycomb_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_honeycomb_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_honeycomb_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_honeycomb_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_honeycomb_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_honeycomb_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_kagome_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_kagome_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_kagome_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_kagome_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_kagome_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_kagome_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_kagome_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_kagome_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_kagome_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_kagome_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_kagome_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_kagome_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_ladder_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_ladder_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_ladder_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_ladder_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_ladder_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_ladder_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_ladder_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_ladder_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_ladder_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_ladder_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_ladder_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_ladder_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_orthorhombic_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_orthorhombic_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_orthorhombic_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_orthorhombic_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_orthorhombic_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_orthorhombic_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_orthorhombic_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_orthorhombic_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_orthorhombic_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_orthorhombic_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_orthorhombic_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_orthorhombic_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_pyrochlore_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_pyrochlore_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_pyrochlore_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_pyrochlore_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_pyrochlore_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_pyrochlore_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_pyrochlore_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_pyrochlore_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_pyrochlore_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_pyrochlore_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_pyrochlore_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_pyrochlore_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_fco_hubbard
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_fco_hubbard ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_fco_hubbardgc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_fco_hubbardgc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_fco_spin
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_fco_spin ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_fco_spingc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_fco_spingc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_fco_kondo
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_fco_kondo ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
add_test(
  NAME matrix_fco_kondogc
  COMMAND ${CMAKE_SOURCE_DIR}/test/uhf/check.sh matrix_fco_kondogc ${CMAKE_SOURCE_DIR}/test/uhf ${TEST_MODE}
)
