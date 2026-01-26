"""
Standard mode for the ladder lattice.

License
-------
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""

from __future__ import annotations

import numpy as np

from stdface_vals import StdIntList
from stdface_model_util import (
    exit_program,
    print_val_d,
    print_val_i,
    init_site,
    set_label,
    not_used_j,
    not_used_c,
    not_used_d,
    not_used_i,
    required_val_i,
    input_spin,
    input_hopp,
    input_coulomb_v,
    malloc_interactions,
    mag_field,
    general_j,
    hubbard_local,
    hopping,
    coulomb,
    print_geometry,
)


def ladder(StdI: StdIntList) -> None:
    """Setup a Hamiltonian for the ladder lattice.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    fp = None
    if StdI.solver != "HWAVE" or StdI.lattice_gp == 1:
        fp = open("lattice.gp", "w")

    # 1. Set lattice size and shape parameters
    print("  @ Lattice Size & Shape\n")

    StdI.a = print_val_d("a", StdI.a, 1.0)
    StdI.length[0] = print_val_d("Wlength", StdI.length[0], StdI.a)
    StdI.length[1] = print_val_d("Llength", StdI.length[1], StdI.a)
    StdI.direct[0, 0] = print_val_d("Wx", StdI.direct[0, 0], StdI.length[0])
    StdI.direct[0, 1] = print_val_d("Wy", StdI.direct[0, 1], 0.0)
    StdI.direct[1, 0] = print_val_d("Lx", StdI.direct[1, 0], 0.0)
    StdI.direct[1, 1] = print_val_d("Ly", StdI.direct[1, 1], StdI.length[1])

    required_val_i("L", StdI.L)
    required_val_i("W", StdI.W)

    not_used_i("a0W", StdI.box[0, 0])
    not_used_i("a0L", StdI.box[0, 1])
    not_used_i("a1W", StdI.box[1, 0])
    not_used_i("a1L", StdI.box[1, 1])

    StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
    not_used_d("phase1", StdI.phase[1])
    StdI.phase[1] = StdI.phase[0]
    StdI.phase[0] = 0.0

    StdI.NsiteUC = StdI.W
    StdI.W = 1
    StdI.direct[0, 0] = float(StdI.NsiteUC)
    init_site(StdI, fp, 2)

    for isite in range(StdI.NsiteUC):
        StdI.tau[isite, 0] = float(isite) / float(StdI.NsiteUC)
        StdI.tau[isite, 1] = 0.0
        StdI.tau[isite, 2] = 0.0

    # 2. Set Hamiltonian parameters
    print("\n  @ Hamiltonian \n")

    not_used_j("J", StdI.JAll, StdI.J)
    not_used_j("J'", StdI.JpAll, StdI.Jp)
    not_used_c("t", StdI.t)
    not_used_c("t'", StdI.tp)
    not_used_d("V", StdI.V)
    not_used_d("V'", StdI.Vp)
    not_used_d("K", StdI.K)

    StdI.h = print_val_d("h", StdI.h, 0.0)
    StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
    StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)

    if StdI.model == "spin":
        StdI.S2 = print_val_i("2S", StdI.S2, 1)
        StdI.D[2, 2] = print_val_d("D", StdI.D[2, 2], 0.0)
        input_spin(StdI.J0, StdI.J0All, "J0")
        input_spin(StdI.J1, StdI.J1All, "J1")
        input_spin(StdI.J2, StdI.J2All, "J2")
        input_spin(StdI.J1p, StdI.J1pAll, "J1'")
        input_spin(StdI.J2p, StdI.J2pAll, "J2'")

        not_used_d("mu", StdI.mu)
        not_used_d("U", StdI.U)
        not_used_c("t0", StdI.t0)
        not_used_c("t1", StdI.t1)
        not_used_c("t2", StdI.t2)
        not_used_c("t1'", StdI.t1p)
        not_used_c("t2'", StdI.t2p)
        not_used_d("V0", StdI.V0)
        not_used_d("V1", StdI.V1)
        not_used_d("V2", StdI.V2)
        not_used_d("V1'", StdI.V1p)
        not_used_d("V2'", StdI.V2p)
    else:
        StdI.mu = print_val_d("mu", StdI.mu, 0.0)
        StdI.U = print_val_d("U", StdI.U, 0.0)
        StdI.t0 = input_hopp(StdI.t, StdI.t0, "t0")
        StdI.t1 = input_hopp(StdI.t, StdI.t1, "t1")
        StdI.t2 = input_hopp(StdI.t, StdI.t2, "t2")
        StdI.t1p = input_hopp(StdI.t, StdI.t1p, "t1'")
        StdI.t2p = input_hopp(StdI.t, StdI.t2p, "t2'")
        StdI.V0 = input_coulomb_v(StdI.V, StdI.V0, "V0")
        StdI.V1 = input_coulomb_v(StdI.V, StdI.V1, "V1")
        StdI.V2 = input_coulomb_v(StdI.V, StdI.V2, "V2")
        StdI.V1p = input_coulomb_v(StdI.V, StdI.V1p, "V1'")
        StdI.V2p = input_coulomb_v(StdI.V, StdI.V2p, "V2'")

        not_used_j("J0", StdI.J0All, StdI.J0)
        not_used_j("J1", StdI.J1All, StdI.J1)
        not_used_j("J2", StdI.J2All, StdI.J2)
        not_used_j("J1p", StdI.J1pAll, StdI.J1p)
        not_used_j("J2p", StdI.J2pAll, StdI.J2p)
        not_used_d("D", StdI.D[2, 2])

        if StdI.model == "hubbard":
            not_used_i("2S", StdI.S2)
            not_used_j("J", StdI.JAll, StdI.J)
        else:
            StdI.S2 = print_val_i("2S", StdI.S2, 1)
            input_spin(StdI.J, StdI.JAll, "J")

    print("\n  @ Numerical conditions\n")

    # 3. Set local spin flags and number of sites
    StdI.nsite = StdI.L * StdI.NsiteUC
    if StdI.model == "kondo":
        StdI.nsite *= 2
    StdI.locspinflag = np.zeros(StdI.nsite, dtype=int)

    if StdI.model == "spin":
        StdI.locspinflag[:] = StdI.S2
    elif StdI.model == "hubbard":
        StdI.locspinflag[:] = 0
    elif StdI.model == "kondo":
        half = StdI.nsite // 2
        StdI.locspinflag[:half] = StdI.S2
        StdI.locspinflag[half:] = 0

    # 4. Calculate maximum number of interactions and allocate arrays
    if StdI.model == "spin":
        ntransMax = StdI.L * StdI.NsiteUC * (StdI.S2 + 1 + 2 * StdI.S2)
        nintrMax = (StdI.L * StdI.NsiteUC * (1 + 1 + 1)
                    * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1)
                    + StdI.L * (StdI.NsiteUC - 1) * (1 + 1 + 1)
                    * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1))
    else:
        ntransMax = (StdI.L * StdI.NsiteUC * 2 * (2 + 2 + 2)
                     + StdI.L * (StdI.NsiteUC - 1) * 2 * (2 + 2 + 2))
        nintrMax = (StdI.L * StdI.NsiteUC * 1
                    + StdI.L * StdI.NsiteUC * 4 * (1 + 1)
                    + StdI.L * (StdI.NsiteUC - 1) * 4 * (1 + 1 + 1))
        if StdI.model == "kondo":
            ntransMax += StdI.L * StdI.NsiteUC * (StdI.S2 + 1 + 2 * StdI.S2)
            nintrMax += StdI.nsite // 2 * (3 * 1 + 1) * (3 * StdI.S2 + 1)

    malloc_interactions(StdI, ntransMax, nintrMax)

    # 5. Set all interactions
    for iL in range(StdI.L):
        for isiteUC in range(StdI.NsiteUC):
            isite = isiteUC + iL * StdI.NsiteUC
            if StdI.model == "kondo":
                isite += StdI.L * StdI.NsiteUC

            # Local terms
            if StdI.model == "spin":
                mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, isite)
                general_j(StdI, StdI.D, StdI.S2, StdI.S2, isite, isite)
            else:
                hubbard_local(StdI, StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, StdI.U, isite)
                if StdI.model == "kondo":
                    jsite = isiteUC + iL * StdI.NsiteUC
                    general_j(StdI, StdI.J, 1, StdI.S2, isite, jsite)
                    mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, jsite)

            # Nearest neighbor along the ladder
            isite, jsite, Cphase, dR = set_label(
                StdI, fp, 0, iL, 0, 1, isiteUC, isiteUC, 1)
            if StdI.model == "spin":
                general_j(StdI, StdI.J1, StdI.S2, StdI.S2, isite, jsite)
            else:
                hopping(StdI, Cphase * StdI.t1, isite, jsite, dR)
                coulomb(StdI, StdI.V1, isite, jsite)

            # Second nearest neighbor along the ladder
            isite, jsite, Cphase, dR = set_label(
                StdI, fp, 0, iL, 0, 2, isiteUC, isiteUC, 2)
            if StdI.model == "spin":
                general_j(StdI, StdI.J1p, StdI.S2, StdI.S2, isite, jsite)
            else:
                hopping(StdI, Cphase * StdI.t1p, isite, jsite, dR)
                coulomb(StdI, StdI.V1p, isite, jsite)

            # Interactions across rungs
            if isiteUC < StdI.NsiteUC - 1:
                # Vertical
                isite, jsite, Cphase, dR = set_label(
                    StdI, fp, 0, iL, 0, 0, isiteUC, isiteUC + 1, 1)
                if StdI.model == "spin":
                    general_j(StdI, StdI.J0, StdI.S2, StdI.S2, isite, jsite)
                else:
                    hopping(StdI, Cphase * StdI.t0, isite, jsite, dR)
                    coulomb(StdI, StdI.V0, isite, jsite)

                # Diagonal 1
                isite, jsite, Cphase, dR = set_label(
                    StdI, fp, 0, iL, 0, 1, isiteUC, isiteUC + 1, 1)
                if StdI.model == "spin":
                    general_j(StdI, StdI.J2, StdI.S2, StdI.S2, isite, jsite)
                else:
                    hopping(StdI, Cphase * StdI.t2, isite, jsite, dR)
                    coulomb(StdI, StdI.V2, isite, jsite)

                # Diagonal 2
                isite, jsite, Cphase, dR = set_label(
                    StdI, fp, 0, iL, 0, -1, isiteUC, isiteUC + 1, 1)
                if StdI.model == "spin":
                    general_j(StdI, StdI.J2p, StdI.S2, StdI.S2, isite, jsite)
                else:
                    hopping(StdI, Cphase * StdI.t2p, isite, jsite, dR)
                    coulomb(StdI, StdI.V2p, isite, jsite)

    if StdI.solver != "HWAVE" or StdI.lattice_gp == 1:
        fp.write("plot '-' w d lc 7\n0.0 0.0\nend\npause -1\n")
        fp.close()

    print_geometry(StdI)


def ladder_boost(StdI: StdIntList) -> None:
    """Setup a boosted Hamiltonian for the ladder lattice (HPhi only).

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    StdI.W = StdI.NsiteUC
    StdI.NsiteUC = 1

    fp = open("boost.def", "w")

    # Magnetic field
    fp.write("# Magnetic field\n")
    fp.write(f"{-0.5 * StdI.Gamma:25.15e} {-0.5 * StdI.Gamma_y:25.15e} {-0.5 * StdI.h:25.15e}\n")

    # Interaction parameters
    fp.write(f"{5}  # Number of type of J\n")

    # J1 - Vertical interactions
    fp.write("# J 1 (inter chain, vertical)\n")
    fp.write(f"{0.25 * StdI.J0[0, 0]:25.15e} {0.25 * StdI.J0[0, 1]:25.15e} {0.25 * StdI.J0[0, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J0[0, 1]:25.15e} {0.25 * StdI.J0[1, 1]:25.15e} {0.25 * StdI.J0[1, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J0[0, 2]:25.15e} {0.25 * StdI.J0[1, 2]:25.15e} {0.25 * StdI.J0[2, 2]:25.15e}\n")

    # J2 - Nearest neighbor along chain
    fp.write("# J 2 (Nearest neighbor, along chain)\n")
    fp.write(f"{0.25 * StdI.J1[0, 0]:25.15e} {0.25 * StdI.J1[0, 1]:25.15e} {0.25 * StdI.J1[0, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J1[0, 1]:25.15e} {0.25 * StdI.J1[1, 1]:25.15e} {0.25 * StdI.J1[1, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J1[0, 2]:25.15e} {0.25 * StdI.J1[1, 2]:25.15e} {0.25 * StdI.J1[2, 2]:25.15e}\n")

    # J3 - Second nearest neighbor along chain
    fp.write("# J 3 (Second nearest neighbor, along chain)\n")
    fp.write(f"{0.25 * StdI.J1p[0, 0]:25.15e} {0.25 * StdI.J1p[0, 1]:25.15e} {0.25 * StdI.J1p[0, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J1p[0, 1]:25.15e} {0.25 * StdI.J1p[1, 1]:25.15e} {0.25 * StdI.J1p[1, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J1p[0, 2]:25.15e} {0.25 * StdI.J1p[1, 2]:25.15e} {0.25 * StdI.J1p[2, 2]:25.15e}\n")

    # J4 - Diagonal 1
    fp.write("# J 4 (inter chain, diagonal1)\n")
    fp.write(f"{0.25 * StdI.J2[0, 0]:25.15e} {0.25 * StdI.J2[0, 1]:25.15e} {0.25 * StdI.J2[0, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J2[0, 1]:25.15e} {0.25 * StdI.J2[1, 1]:25.15e} {0.25 * StdI.J2[1, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J2[0, 2]:25.15e} {0.25 * StdI.J2[1, 2]:25.15e} {0.25 * StdI.J2[2, 2]:25.15e}\n")

    # J5 - Diagonal 2
    fp.write("# J 5 (inter chain, diagonal2)\n")
    fp.write(f"{0.25 * StdI.J2p[0, 0]:25.15e} {0.25 * StdI.J2p[0, 1]:25.15e} {0.25 * StdI.J2p[0, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J2p[0, 1]:25.15e} {0.25 * StdI.J2p[1, 1]:25.15e} {0.25 * StdI.J2p[1, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J2p[0, 2]:25.15e} {0.25 * StdI.J2p[1, 2]:25.15e} {0.25 * StdI.J2p[2, 2]:25.15e}\n")

    # Validate parameters
    if StdI.S2 != 1:
        print("\n ERROR! S2 must be 1 in Boost. \n")
        exit_program(-1)
    StdI.ishift_nspin = 2
    if StdI.W != 2:
        print("\n ERROR! W != 2 \n")
        exit_program(-1)
    if StdI.L % 2 != 0:
        print("\n ERROR! L %% 2 != 0 \n")
        exit_program(-1)
    if StdI.L < 4:
        print("\n ERROR! L < 4 \n")
        exit_program(-1)

    StdI.W = StdI.L
    StdI.L = 2
    StdI.num_pivot = StdI.W // 2

    fp.write("# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n")
    fp.write(f"{StdI.W} {StdI.L} {StdI.num_pivot} {StdI.ishift_nspin}\n")

    # 6-spin star list
    StdI.list_6spin_star = np.zeros((StdI.num_pivot, 7), dtype=int)
    for ipivot in range(StdI.num_pivot):
        StdI.list_6spin_star[ipivot, 0] = 7  # num of J
        StdI.list_6spin_star[ipivot, 1] = 1
        StdI.list_6spin_star[ipivot, 2] = 1
        StdI.list_6spin_star[ipivot, 3] = 1
        StdI.list_6spin_star[ipivot, 4] = 1
        StdI.list_6spin_star[ipivot, 5] = 1
        StdI.list_6spin_star[ipivot, 6] = 1  # flag

    fp.write("# StdI->list_6spin_star\n")
    for ipivot in range(StdI.num_pivot):
        fp.write(f"# pivot {ipivot}\n")
        for isite in range(7):
            fp.write(f"{StdI.list_6spin_star[ipivot, isite]} ")
        fp.write("\n")

    # 6-spin pair list
    StdI.list_6spin_pair = np.zeros((StdI.num_pivot, 7, 7), dtype=int)
    for ipivot in range(StdI.num_pivot):
        # kintr=0
        StdI.list_6spin_pair[ipivot, :, 0] = [0, 1, 2, 3, 4, 5, 1]
        # kintr=1
        StdI.list_6spin_pair[ipivot, :, 1] = [0, 2, 1, 3, 4, 5, 2]
        # kintr=2
        StdI.list_6spin_pair[ipivot, :, 2] = [1, 3, 0, 2, 4, 5, 2]
        # kintr=3
        StdI.list_6spin_pair[ipivot, :, 3] = [0, 4, 1, 2, 3, 5, 3]
        # kintr=4
        StdI.list_6spin_pair[ipivot, :, 4] = [1, 5, 0, 2, 3, 4, 3]
        # kintr=5
        StdI.list_6spin_pair[ipivot, :, 5] = [0, 3, 1, 2, 4, 5, 4]
        # kintr=6
        StdI.list_6spin_pair[ipivot, :, 6] = [1, 2, 0, 3, 4, 5, 5]

    fp.write("# StdI->list_6spin_pair\n")
    for ipivot in range(StdI.num_pivot):
        fp.write(f"# pivot {ipivot}\n")
        for kintr in range(StdI.list_6spin_star[ipivot, 0]):
            for isite in range(7):
                fp.write(f"{StdI.list_6spin_pair[ipivot, isite, kintr]} ")
            fp.write("\n")

    fp.close()
