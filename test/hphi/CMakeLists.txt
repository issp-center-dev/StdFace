add_test(
  NAME ctpq_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh ctpq_spin_kagome ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME ctpq_spin_kagome_randomsphere
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh ctpq_spin_kagome_randomsphere ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_genspin_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_genspin_ladder ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_genspingc_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_genspingc_ladder ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_ham_io
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_ham_io ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_hubbard_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_hubbard_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_hubbardgc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_hubbardgc_tri ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_kondo_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_kondo_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_kondogc_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_kondogc_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_spin_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_spin_tri ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME fulldiag_spingc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh fulldiag_spingc_tri ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_genspin_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_genspin_ladder ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_genspin_ladder_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_genspin_ladder_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_genspingc_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_genspingc_ladder ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_hubbard_square ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_hubbard_square_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_hubbard_square_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_hubbardgc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_hubbardgc_tri ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_kondo_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_kondo_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_kondo_chain_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_kondo_chain_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_kondogc_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_kondogc_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_spinGC_DM
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_spinGC_DM ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_spin_kagome ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_spin_kagome_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_spin_kagome_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_spingc_hcor
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_spingc_hcor ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lanczos_spingc_honey
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lanczos_spingc_honey ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_genspin_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_genspin_ladder ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_genspin_ladder_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_genspin_ladder_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_genspingc_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_genspingc_ladder ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_hubbard_square ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_hubbard_square_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_hubbard_square_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_hubbard_szconserv
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_hubbard_szconserv ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_hubbardgc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_hubbardgc_tri ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_kondo_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_kondo_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_kondo_chain_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_kondo_chain_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_kondogc_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_kondogc_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_spin_kagome ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_spin_kagome_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_spin_kagome_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_spingc_6body
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_spingc_6body ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_spingc_Sy
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_spingc_Sy ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_spingc_hcor
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_spingc_hcor ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME lobcg_spingc_honey
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh lobcg_spingc_honey ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME spectrum_genspin_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh spectrum_genspin_ladder ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME spectrum_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh spectrum_hubbard_square ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME spectrum_hubbardgc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh spectrum_hubbardgc_tri ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME spectrum_kondo_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh spectrum_kondo_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME spectrum_kondogc_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh spectrum_kondogc_chain ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME spectrum_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh spectrum_spin_kagome ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME spectrum_spingc_honey
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh spectrum_spingc_honey ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_ac_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_ac_hubbard_square ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_ac_hubbard_square_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_ac_hubbard_square_restart ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_dc_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_dc_hubbard_square ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_hubbard_chain_interall
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_hubbard_chain_interall ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_hubbard_chain_interall_diagonal
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_hubbard_chain_interall_diagonal ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_kondo_chain_interall
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_kondo_chain_interall ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_pulse_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_pulse_hubbard_square ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_quench_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_quench_hubbard_square ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME te_spin_chain_interall
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh te_spin_chain_interall ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME tpq_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh tpq_spin_kagome ${CMAKE_SOURCE_DIR}/test/hphi
)
add_test(
  NAME tpq_spin_kagome_randomsphere
  COMMAND ${CMAKE_SOURCE_DIR}/test/hphi/check.sh tpq_spin_kagome_randomsphere ${CMAKE_SOURCE_DIR}/test/hphi
)
