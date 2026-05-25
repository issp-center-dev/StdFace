set(TEST_MODE "base")
#set(TEST_MODE "python")

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
