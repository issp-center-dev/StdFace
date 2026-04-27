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

# --- matrix tests (auto-generated) ---
add_test(
  NAME matrix_chain_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_chain_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_chain_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_square_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_square_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_triangular_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_triangular_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_honeycomb_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_honeycomb_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_kagome_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_kagome_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_ladder_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_ladder_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_orthorhombic_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_orthorhombic_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_pyrochlore_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_pyrochlore_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_kondo_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_kondo_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME matrix_fco_kondogc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh matrix_fco_kondogc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)

# --- wannier90 tests (auto-generated) ---
add_test(
  NAME wannier90_hubbard_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh wannier90_hubbard_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME wannier90_hubbard_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh wannier90_hubbard_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME wannier90_hubbardgc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh wannier90_hubbardgc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME wannier90_hubbardgc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh wannier90_hubbardgc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME wannier90_spin_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh wannier90_spin_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME wannier90_spin_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh wannier90_spin_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME wannier90_spingc_uhfr
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh wannier90_spingc_uhfr ${CMAKE_SOURCE_DIR}/test/hwave
)
add_test(
  NAME wannier90_spingc_uhfk
  COMMAND ${CMAKE_SOURCE_DIR}/test/hwave/check.sh wannier90_spingc_uhfk ${CMAKE_SOURCE_DIR}/test/hwave
)
