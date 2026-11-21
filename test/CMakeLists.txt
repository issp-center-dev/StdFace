set(TEST_MODE "python" CACHE STRING "Integration test mode: python (stdface CLI) or base (C dry binaries)")

if (HPHI)
    add_subdirectory(hphi)
endif()

if (MVMC)
    add_subdirectory(mvmc)
endif()

if (UHF)
    add_subdirectory(uhf)
endif()

if (HWAVE)
    add_subdirectory(hwave)
endif()

# The tests below exercise the C generators directly (linked library or
# dry-run binary), so they are registered for the base mode only.
if(TEST_MODE STREQUAL "base")
foreach(mode hphi mvmc uhf hwave)
    if(TARGET StdFace_${mode})
        add_executable(transfer_capacity_${mode} transfer_capacity.c)
        target_include_directories(transfer_capacity_${mode} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
        target_link_libraries(transfer_capacity_${mode} StdFace_${mode} m)
        add_test(NAME TransferCapacity_${mode} COMMAND transfer_capacity_${mode})
        foreach(invalid negative-transfers negative-interactions negative-sites overflow)
            if(invalid MATCHES "^negative-(transfers|interactions)$")
                set(expected "Error: negative interaction allocation size.")
            else()
                set(expected "Error: local transfer allocation size is invalid or too large.")
            endif()
            add_test(NAME TransferCapacity_${mode}_${invalid}
                COMMAND ${CMAKE_COMMAND}
                    -DPROGRAM=$<TARGET_FILE:transfer_capacity_${mode}>
                    -DARGUMENT=${invalid}
                    "-DEXPECTED=${expected}"
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_transfer_capacity_failure.cmake)
        endforeach()
    endif()
endforeach()

foreach(mode hphi mvmc uhf hwave)
    if(TARGET ${mode}_dry.out)
        add_test(NAME Unused3DParameters_${mode}
            COMMAND ${CMAKE_COMMAND}
                -DPROGRAM=$<TARGET_FILE:${mode}_dry.out>
                -DMODE=${mode}
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/unused-3d-${mode}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/check_unused_3d.cmake)
        set_tests_properties(Unused3DParameters_${mode} PROPERTIES TIMEOUT 180)
        add_test(NAME WannierSpinU_${mode}
            COMMAND ${CMAKE_COMMAND}
                -DPROGRAM=$<TARGET_FILE:${mode}_dry.out>
                -DMODE=${mode}
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/wannier-spin-u-${mode}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/check_wannier_spin_u.cmake)
        set_tests_properties(WannierSpinU_${mode} PROPERTIES TIMEOUT 180)
        add_test(NAME WannierInitial_${mode}
            COMMAND ${CMAKE_COMMAND}
                -DPROGRAM=$<TARGET_FILE:${mode}_dry.out>
                -DMODE=${mode}
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/wannier-initial-${mode}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/check_wannier_initial.cmake)
        set_tests_properties(WannierInitial_${mode} PROPERTIES TIMEOUT 180)
        add_test(NAME WannierReader_${mode}
            COMMAND ${CMAKE_COMMAND}
                -DPROGRAM=$<TARGET_FILE:${mode}_dry.out>
                -DMODE=${mode}
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/wannier-reader-${mode}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/check_wannier_reader.cmake)
        set_tests_properties(WannierReader_${mode} PROPERTIES TIMEOUT 180)
        add_test(NAME WannierDensity_${mode}
            COMMAND ${CMAKE_COMMAND}
                -DPROGRAM=$<TARGET_FILE:${mode}_dry.out>
                -DMODE=${mode}
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/wannier-density-${mode}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/check_wannier_density.cmake)
        set_tests_properties(WannierDensity_${mode} PROPERTIES TIMEOUT 180)
    endif()
endforeach()
endif()

# In python mode, the generator scripts above run against the Python CLI
# through an exit-status adapter. check_wannier_reader.cmake is left out: the
# Python reader does not mirror the C reader's full structural validation.
if(TEST_MODE STREQUAL "python")
foreach(pair hphi:HPhi mvmc:mVMC uhf:UHF hwave:HWAVE)
    string(REPLACE ":" ";" pair "${pair}")
    list(GET pair 0 mode)
    list(GET pair 1 solver)
    string(TOUPPER ${mode} option)
    if(${option})
        foreach(check Unused3DParameters:unused_3d WannierSpinU:wannier_spin_u
                WannierInitial:wannier_initial WannierDensity:wannier_density)
            string(REPLACE ":" ";" check "${check}")
            list(GET check 0 name)
            list(GET check 1 script)
            add_test(NAME ${name}_${mode}
                COMMAND ${CMAKE_COMMAND} -E env STDFACE_SOLVER=${solver}
                    ${CMAKE_COMMAND}
                    -DPROGRAM=${CMAKE_CURRENT_SOURCE_DIR}/common/python_cli.sh
                    -DMODE=${mode}
                    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/${script}-${mode}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_${script}.cmake)
            set_tests_properties(${name}_${mode} PROPERTIES TIMEOUT 900)
        endforeach()
    endif()
endforeach()
endif()
