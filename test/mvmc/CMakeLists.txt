add_test(
  NAME HeisenbergChain
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HeisenbergChain ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME HeisenbergChain_cmp
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HeisenbergChain_cmp ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME HeisenbergChain_fsz
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HeisenbergChain_fsz ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME HubbardChain
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HubbardChain ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME HubbardChainLanczos
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HubbardChainLanczos ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME HubbardChain_cmp
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HubbardChain_cmp ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME HubbardChain_fsz
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HubbardChain_fsz ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME HubbardTetragonal
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HubbardTetragonal ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME HubbardTetragonal_MomentumProjection
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh HubbardTetragonal_MomentumProjection ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME KondoChain
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh KondoChain ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME KondoChain_Stot1_cmp
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh KondoChain_Stot1_cmp ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME KondoChain_cmp
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh KondoChain_cmp ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME KondoChain_fsz
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh KondoChain_fsz ${CMAKE_SOURCE_DIR}/test/mvmc
)
add_test(
  NAME SpinChainLanczos
  COMMAND ${CMAKE_SOURCE_DIR}/test/mvmc/check.sh SpinChainLanczos ${CMAKE_SOURCE_DIR}/test/mvmc
)
