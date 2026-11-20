include("${CMAKE_CURRENT_LIST_DIR}/wannier_fixture.cmake")

# A density file need not list every R of the interaction files: an unlisted R
# has zero density. Each input below must therefore generate the same
# definition files as the input whose density file lists that R as zero.
function(density_input dir n dc)
    wannier_input("${dir}" ${n} Hubbard
        "doublecounting = \"${dc}\"\ncutoff_length_U = -1\ncutoff_length_J = -1\n")
    file(READ "${dir}/stan.in" input)
    string(REPLACE "L = 1\nHeight = 1\n" "L = 3\nHeight = 3\n" cubic "${input}")
    if(cubic STREQUAL input)
        message(FATAL_ERROR "Invalid test fixture ${dir}: lattice size not found")
    endif()
    file(WRITE "${dir}/stan.in" "${cubic}")
endfunction()

# Compare the values of the generated rows, as "real|imaginary", with the
# expected ones; with offdiagonal, only rows between two different sites.
function(expect_values label path offdiagonal expected)
    file(STRINGS "${path}" rows REGEX "^[ ]*[0-9]+[ ]+[0-9]+")
    foreach(row ${rows})
        string(STRIP "${row}" row)
        string(REPLACE "-0.000000000000000" "0.000000000000000" row "${row}")
        string(REGEX REPLACE " +" ";" row "${row}")
        list(GET row 0 i)
        list(GET row 2 j)
        if(NOT offdiagonal OR NOT i EQUAL j)
            list(GET row 4 real)
            list(GET row 5 imaginary)
            list(FIND expected "${real}|${imaginary}" index)
            if(index EQUAL -1)
                message(FATAL_ERROR "${label}: unexpected entry ${row} in ${path}")
            endif()
            list(REMOVE_AT expected ${index})
        endif()
    endforeach()
    if(expected)
        message(FATAL_ERROR "${label}: missing entries in ${path}: ${expected}")
    endif()
endfunction()

# Each value once per cell and spin: 27 cells, two spins.
function(per_cell_and_spin out)
    set(values "")
    foreach(value ${ARGN})
        foreach(k RANGE 1 54)
            list(APPEND values "${value}")
        endforeach()
    endforeach()
    set(${out} "${values}" PARENT_SCOPE)
endfunction()

# One orbital, a single nonlocal U or J term at R, density at R = 0 only: the
# lookup is below or above the density file on each axis.
set(count 0)
per_cell_and_spin(onsite "0.250000000000000|0.000000000000000")
foreach(family u j)
    foreach(axis 0 1 2)
        foreach(side below above)
            set(r "")
            foreach(k 0 1 2)
                if(NOT k EQUAL axis)
                    list(APPEND r 0)
                elseif(side STREQUAL "below")
                    list(APPEND r -1)
                else()
                    list(APPEND r 1)
                endif()
            endforeach()
            string(REPLACE ";" " " r "${r}")
            set(base "${WORK}/${family}-${axis}-${side}")
            foreach(dc hartree hartree_u full)
                foreach(density unlisted zero)
                    set(dir "${base}-${dc}-${density}")
                    density_input("${dir}" 1 ${dc})
                    write_wannier("${dir}/zvo_hr.dat" 1 "-1 0 0;0 0 0;1 0 0" "-0.5 0;-0.2 0;-0.5 0")
                    if(family STREQUAL "u")
                        write_wannier("${dir}/zvo_ur.dat" 1 "0 0 0;${r}" "4 0;0.25 0")
                        write_wannier("${dir}/zvo_jr.dat" 1 "0 0 0" "0 0")
                    else()
                        write_wannier("${dir}/zvo_ur.dat" 1 "0 0 0" "4 0")
                        write_wannier("${dir}/zvo_jr.dat" 1 "0 0 0;${r}" "0 0;0.25 0")
                    endif()
                    if(density STREQUAL "unlisted")
                        write_wannier("${dir}/zvo_dr.dat" 1 "0 0 0" "0.5 0")
                    else()
                        write_wannier("${dir}/zvo_dr.dat" 1 "0 0 0;${r}" "0.5 0;0 0")
                    endif()
                    run_wannier("${dir}" FALSE "")
                    math(EXPR count "${count} + 1")
                endforeach()
                same_def("${base}-${dc}-zero" "${base}-${dc}-unlisted")
                # The initial guess has the on-site density only (half of it).
                expect_values("${family}-${axis}-${side}-${dc}"
                    "${base}-${dc}-unlisted/initial.def" FALSE "${onsite}")
            endforeach()
            # Without density at R the exchange (Fock) correction vanishes.
            same_def("${base}-hartree-unlisted" "${base}-full-unlisted")
        endforeach()
    endforeach()
endforeach()

# Two orbitals with listed and unlisted R in one input.  The density file spans
# R0 = 0..1 only.  Two U terms are above it on axis 1 and below it on axis 2,
# at offsets that would alias a listed block if a bound were ignored.  The U and
# J terms at R = 1 0 0 use different, listed off-diagonal elements.
foreach(dc hartree full)
    foreach(density unlisted zero)
        set(dir "${WORK}/mixed-${dc}-${density}")
        density_input("${dir}" 2 ${dc})
        write_wannier("${dir}/zvo_hr.dat" 2 "0 0 0" "-0.2 0;0 0;0 0;-0.2 0")
        write_wannier("${dir}/zvo_ur.dat" 2 "0 0 0;1 0 0;0 1 0;1 0 -1"
            "4 0;0 0;0 0;8 0;0 0;0.25 0;0 0;0 0;0 0;0 0;0.3 0;0 0;0.35 0;0 0;0 0;0 0")
        write_wannier("${dir}/zvo_jr.dat" 2 "0 0 0;1 0 0" "0 0;0 0;0 0;0 0;0 0;0 0;0.2 0;0 0")
        set(rs "0 0 0;1 0 0")
        set(values "0.5 0;0 0;0 0;0.75 0;0.11 0.01;0.12 0.02;0.13 0.03;0.14 0.04")
        if(density STREQUAL "zero")
            set(rs "${rs};0 1 0;1 0 -1")
            set(values "${values};0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0")
        endif()
        write_wannier("${dir}/zvo_dr.dat" 2 "${rs}" "${values}")
        run_wannier("${dir}" FALSE "")
        math(EXPR count "${count} + 1")
    endforeach()
    same_def("${WORK}/mixed-${dc}-zero" "${WORK}/mixed-${dc}-unlisted")
    # Half the on-site densities, of D(1 0 0)[1][2] = 0.12 + 0.02i for the U
    # term and of D(1 0 0)[2][1] = 0.13 + 0.03i for the J term, each with its
    # conjugate for the reverse pair. The unlisted elements add no entry.
    per_cell_and_spin(expected
        "0.250000000000000|0.000000000000000" "0.375000000000000|0.000000000000000"
        "0.060000000000000|0.010000000000000" "0.060000000000000|-0.010000000000000"
        "0.065000000000000|0.015000000000000" "0.065000000000000|-0.015000000000000")
    expect_values("mixed-${dc}" "${WORK}/mixed-${dc}-unlisted/initial.def" FALSE "${expected}")
endforeach()
# The only transfers between sites are the exchange (Fock) corrections of the
# listed elements: -U D / 2 with U = 0.25, and J (D + 2 Re D) / 2 with J = 0.2.
expect_values("mixed-hartree" "${WORK}/mixed-hartree-unlisted/trans.def" TRUE "")
per_cell_and_spin(expected
    "-0.015000000000000|-0.002500000000000" "-0.015000000000000|0.002500000000000"
    "0.039000000000000|0.003000000000000" "0.039000000000000|-0.003000000000000")
expect_values("mixed-full" "${WORK}/mixed-full-unlisted/trans.def" TRUE "${expected}")
message(STATUS "Wannier density domain (${MODE}): ${count} inputs passed")
