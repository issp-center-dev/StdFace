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

foreach(mode hphi mvmc uhf hwave)
    if(TARGET ${mode}_dry.out)
        add_test(NAME Unused3DParameters_${mode}
            COMMAND ${CMAKE_COMMAND}
                -DPROGRAM=$<TARGET_FILE:${mode}_dry.out>
                -DMODE=${mode}
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/unused-3d-${mode}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/check_unused_3d.cmake)
        set_tests_properties(Unused3DParameters_${mode} PROPERTIES TIMEOUT 180)
    endif()
endforeach()
