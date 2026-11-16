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
    if(TARGET StdFace_${mode})
        add_executable(transfer_capacity_${mode} transfer_capacity.c)
        target_include_directories(transfer_capacity_${mode} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
        target_link_libraries(transfer_capacity_${mode} StdFace_${mode} m)
        add_test(NAME TransferCapacity_${mode} COMMAND transfer_capacity_${mode})
        foreach(invalid negative-transfers negative-interactions negative-sites overflow)
            add_test(NAME TransferCapacity_${mode}_${invalid}
                COMMAND transfer_capacity_${mode} ${invalid})
            set_tests_properties(TransferCapacity_${mode}_${invalid} PROPERTIES WILL_FAIL TRUE)
        endforeach()
    endif()
endforeach()
