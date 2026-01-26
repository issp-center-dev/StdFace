"""
Standard mode for the kagome lattice.

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

import math

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
    input_spin_nn,
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


def kagome(StdI: StdIntList) -> None:
    """Setup a Hamiltonian for the kagome lattice.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    # (1) Compute the shape of the super-cell and sites in the super-cell
    fp = None
    if StdI.solver != "HWAVE" or StdI.lattice_gp == 1:
        fp = open("lattice.gp", "w")

    StdI.NsiteUC = 3

    print("  @ Lattice Size & Shape\n")

    StdI.a = print_val_d("a", StdI.a, 1.0)
    StdI.length[0] = print_val_d("Wlength", StdI.length[0], StdI.a)
    StdI.length[1] = print_val_d("Llength", StdI.length[1], StdI.a)
    StdI.direct[0, 0] = print_val_d("Wx", StdI.direct[0, 0], StdI.length[0])
    StdI.direct[0, 1] = print_val_d("Wy", StdI.direct[0, 1], 0.0)
    StdI.direct[1, 0] = print_val_d("Lx", StdI.direct[1, 0], StdI.length[1] * 0.5)
    StdI.direct[1, 1] = print_val_d("Ly", StdI.direct[1, 1], StdI.length[1] * 0.5 * math.sqrt(3.0))

    StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
    StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)

    init_site(StdI, fp, 2)
    StdI.tau[0, 0] = 0.0
    StdI.tau[0, 1] = 0.0
    StdI.tau[0, 2] = 0.0
    StdI.tau[1, 0] = 0.5
    StdI.tau[1, 1] = 0.0
    StdI.tau[1, 2] = 0.0
    StdI.tau[2, 0] = 0.0
    StdI.tau[2, 1] = 0.5
    StdI.tau[2, 2] = 0.0

    # (2) check & store parameters of Hamiltonian
    print("\n  @ Hamiltonian \n")
    not_used_d("K", StdI.K)
    StdI.h = print_val_d("h", StdI.h, 0.0)
    StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
    StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)

    if StdI.model == "spin":
        StdI.S2 = print_val_i("2S", StdI.S2, 1)
        StdI.D[2, 2] = print_val_d("D", StdI.D[2, 2], 0.0)
        input_spin_nn(StdI.J, StdI.JAll, StdI.J0, StdI.J0All, "J0")
        input_spin_nn(StdI.J, StdI.JAll, StdI.J1, StdI.J1All, "J1")
        input_spin_nn(StdI.J, StdI.JAll, StdI.J2, StdI.J2All, "J2")
        input_spin_nn(StdI.Jp, StdI.JpAll, StdI.J0p, StdI.J0pAll, "J0'")
        input_spin_nn(StdI.Jp, StdI.JpAll, StdI.J1p, StdI.J1pAll, "J1'")
        input_spin_nn(StdI.Jp, StdI.JpAll, StdI.J2p, StdI.J2pAll, "J2'")

        not_used_d("mu", StdI.mu)
        not_used_d("U", StdI.U)
        not_used_c("t", StdI.t)
        not_used_c("t0", StdI.t)
        not_used_c("t0", StdI.t0)
        not_used_c("t1", StdI.t1)
        not_used_c("t2", StdI.t2)
        not_used_c("t'", StdI.tp)
        not_used_c("t0'", StdI.t0p)
        not_used_c("t1'", StdI.t1p)
        not_used_c("t2'", StdI.t2p)
        not_used_d("V", StdI.V)
        not_used_d("V0", StdI.V0)
        not_used_d("V1", StdI.V1)
        not_used_d("V2", StdI.V2)
        not_used_d("V'", StdI.Vp)
        not_used_d("V0'", StdI.V0p)
        not_used_d("V1'", StdI.V1p)
        not_used_d("V2'", StdI.V2p)
    else:
        StdI.mu = print_val_d("mu", StdI.mu, 0.0)
        StdI.U = print_val_d("U", StdI.U, 0.0)
        StdI.t0 = input_hopp(StdI.t, StdI.t0, "t0")
        StdI.t1 = input_hopp(StdI.t, StdI.t1, "t1")
        StdI.t2 = input_hopp(StdI.t, StdI.t2, "t2")
        StdI.t0p = input_hopp(StdI.tp, StdI.t0p, "t0'")
        StdI.t1p = input_hopp(StdI.tp, StdI.t1p, "t1'")
        StdI.t2p = input_hopp(StdI.tp, StdI.t2p, "t2'")
        StdI.V0 = input_coulomb_v(StdI.V, StdI.V0, "V0")
        StdI.V1 = input_coulomb_v(StdI.V, StdI.V1, "V1")
        StdI.V2 = input_coulomb_v(StdI.V, StdI.V2, "V2")
        StdI.V0p = input_coulomb_v(StdI.Vp, StdI.V0p, "V0'")
        StdI.V1p = input_coulomb_v(StdI.Vp, StdI.V1p, "V1'")
        StdI.V2p = input_coulomb_v(StdI.Vp, StdI.V2p, "V2'")

        not_used_j("J0", StdI.J0All, StdI.J0)
        not_used_j("J1", StdI.J1All, StdI.J1)
        not_used_j("J2", StdI.J2All, StdI.J2)
        not_used_j("J'", StdI.JpAll, StdI.Jp)
        not_used_j("J0'", StdI.J0pAll, StdI.J0p)
        not_used_j("J1'", StdI.J1pAll, StdI.J1p)
        not_used_j("J2'", StdI.J2pAll, StdI.J2p)
        not_used_d("D", StdI.D[2, 2])

        if StdI.model == "hubbard":
            not_used_i("2S", StdI.S2)
            not_used_j("J", StdI.JAll, StdI.J)
        else:
            StdI.S2 = print_val_i("2S", StdI.S2, 1)
            input_spin(StdI.J, StdI.JAll, "J")

    print("\n  @ Numerical conditions\n")

    # (3) Set local spin flag and number of sites
    StdI.nsite = StdI.NsiteUC * StdI.NCell
    if StdI.model == "kondo":
        StdI.nsite *= 2
    StdI.locspinflag = np.zeros(StdI.nsite, dtype=int)

    if StdI.model == "spin":
        StdI.locspinflag[:] = StdI.S2
    elif StdI.model == "hubbard":
        StdI.locspinflag[:] = 0
    else:
        half = StdI.nsite // 2
        StdI.locspinflag[:half] = StdI.S2
        StdI.locspinflag[half:] = 0

    # (4) Compute upper limit of Transfer & Interaction
    if StdI.model == "spin":
        ntransMax = StdI.nsite * (StdI.S2 + 1 + 2 * StdI.S2)
        nintrMax = (StdI.NCell * (StdI.NsiteUC + 6 + 6)
                    * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1))
    else:
        ntransMax = StdI.NCell * 2 * (2 * StdI.NsiteUC + 12 + 12)
        nintrMax = StdI.NCell * (StdI.NsiteUC + 4 * (6 + 6))
        if StdI.model == "kondo":
            ntransMax += StdI.nsite // 2 * (StdI.S2 + 1 + 2 * StdI.S2)
            nintrMax += StdI.nsite // 2 * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1)

    malloc_interactions(StdI, ntransMax, nintrMax)

    # (5) Set Transfer & Interaction
    for kCell in range(StdI.NCell):
        iW = StdI.Cell[kCell, 0]
        iL = StdI.Cell[kCell, 1]

        # Local term
        isite = StdI.NsiteUC * kCell
        if StdI.model == "kondo":
            isite += StdI.nsite // 2

        if StdI.model == "spin":
            for isiteUC in range(StdI.NsiteUC):
                mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, isite + isiteUC)
                general_j(StdI, StdI.D, StdI.S2, StdI.S2, isite + isiteUC, isite + isiteUC)
        else:
            for isiteUC in range(StdI.NsiteUC):
                hubbard_local(StdI, StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, StdI.U, isite + isiteUC)
            if StdI.model == "kondo":
                jsite = StdI.NsiteUC * kCell
                for isiteUC in range(StdI.NsiteUC):
                    general_j(StdI, StdI.J, 1, StdI.S2, isite + isiteUC, jsite + isiteUC)
                    mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, jsite + isiteUC)

        # Nearest neighbor intra cell 0 -> 1
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 0, 0, 0, 1, 1)
        if StdI.model == "spin":
            general_j(StdI, StdI.J2, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t2, isite, jsite, dR)
            coulomb(StdI, StdI.V2, isite, jsite)

        # Nearest neighbor intra cell 0 -> 2
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 0, 0, 0, 2, 1)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1, isite, jsite, dR)
            coulomb(StdI, StdI.V1, isite, jsite)

        # Nearest neighbor intra cell 1 -> 2
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 0, 0, 1, 2, 1)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0, isite, jsite, dR)
            coulomb(StdI, StdI.V0, isite, jsite)

        # Nearest neighbor along W 1 -> 0
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 1, 0, 1, 0, 1)
        if StdI.model == "spin":
            general_j(StdI, StdI.J2, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t2, isite, jsite, dR)
            coulomb(StdI, StdI.V2, isite, jsite)

        # Nearest neighbor along L 2 -> 0
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 0, 1, 2, 0, 1)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1, isite, jsite, dR)
            coulomb(StdI, StdI.V1, isite, jsite)

        # Nearest neighbor along W-L 1 -> 2
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 1, -1, 1, 2, 1)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0, isite, jsite, dR)
            coulomb(StdI, StdI.V0, isite, jsite)

        # Second nearest neighbor along W 2 -> 0
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 1, 0, 2, 0, 2)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1p, isite, jsite, dR)
            coulomb(StdI, StdI.V1p, isite, jsite)

        # Second nearest neighbor along W 1 -> 2
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 1, 0, 1, 2, 2)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0p, isite, jsite, dR)
            coulomb(StdI, StdI.V0p, isite, jsite)

        # Second nearest neighbor along L 1 -> 0
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 0, 1, 1, 0, 2)
        if StdI.model == "spin":
            general_j(StdI, StdI.J2p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t2p, isite, jsite, dR)
            coulomb(StdI, StdI.V2p, isite, jsite)

        # Second nearest neighbor along L 2 -> 1
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 0, 1, 2, 1, 2)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0p, isite, jsite, dR)
            coulomb(StdI, StdI.V0p, isite, jsite)

        # Second nearest neighbor along W-L 0 -> 2
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 1, -1, 0, 2, 2)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1p, isite, jsite, dR)
            coulomb(StdI, StdI.V1p, isite, jsite)

        # Second nearest neighbor along L-W 0 -> 1
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, -1, 1, 0, 1, 2)
        if StdI.model == "spin":
            general_j(StdI, StdI.J2p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t2p, isite, jsite, dR)
            coulomb(StdI, StdI.V2p, isite, jsite)

    if StdI.solver != "HWAVE" or StdI.lattice_gp == 1:
        fp.write("plot '-' w d lc 7\n0.0 0.0\nend\npause -1\n")
        fp.close()

    print_geometry(StdI)


def kagome_boost(StdI: StdIntList) -> None:
    """Setup a boosted Hamiltonian for the kagome lattice (HPhi only).

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    if StdI.box[0, 1] != 0 or StdI.box[1, 0] != 0:
        print("\nERROR ! (a0W, a0L, a1W, a1L) can not be used with SpinGCBoost.\n")
        exit_program(-1)

    for i1 in range(3):
        for i2 in range(3):
            if abs(StdI.Jp[i1, i2]) > 1.0e-8:
                print("\nERROR ! J' can not be used with SpinGCBoost.\n")
                exit_program(-1)

    # Magnetic field
    fp = open("boost.def", "w")
    fp.write("# Magnetic field\n")
    fp.write(f"{-0.5 * StdI.Gamma:25.15e} {-0.5 * StdI.Gamma_y:25.15e} {-0.5 * StdI.h:25.15e}\n")

    # Interaction
    fp.write(f"{3}  # Number of type of J\n")
    fp.write("# J 0\n")
    fp.write(f"{0.25 * StdI.J0[0, 0]:25.15e} {0.25 * StdI.J0[0, 1]:25.15e} {0.25 * StdI.J0[0, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J0[0, 1]:25.15e} {0.25 * StdI.J0[1, 1]:25.15e} {0.25 * StdI.J0[1, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J0[0, 2]:25.15e} {0.25 * StdI.J0[1, 2]:25.15e} {0.25 * StdI.J0[2, 2]:25.15e}\n")
    fp.write("# J 1\n")
    fp.write(f"{0.25 * StdI.J1[0, 0]:25.15e} {0.25 * StdI.J1[0, 1]:25.15e} {0.25 * StdI.J1[0, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J1[0, 1]:25.15e} {0.25 * StdI.J1[1, 1]:25.15e} {0.25 * StdI.J1[1, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J1[0, 2]:25.15e} {0.25 * StdI.J1[1, 2]:25.15e} {0.25 * StdI.J1[2, 2]:25.15e}\n")
    fp.write("# J 2\n")
    fp.write(f"{0.25 * StdI.J2[0, 0]:25.15e} {0.25 * StdI.J2[0, 1]:25.15e} {0.25 * StdI.J2[0, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J2[0, 1]:25.15e} {0.25 * StdI.J2[1, 1]:25.15e} {0.25 * StdI.J2[1, 2]:25.15e}\n")
    fp.write(f"{0.25 * StdI.J2[0, 2]:25.15e} {0.25 * StdI.J2[1, 2]:25.15e} {0.25 * StdI.J2[2, 2]:25.15e}\n")

    # Topology
    if StdI.S2 != 1:
        print("\n ERROR! S2 must be 1 in Boost. \n")
        exit_program(-1)
    StdI.ishift_nspin = 3
    if StdI.L < 2:
        print("\n ERROR! L < 2 \n")
        exit_program(-1)
    if StdI.W % StdI.ishift_nspin != 0:
        print(f"\n ERROR! W %% {StdI.ishift_nspin} != 0 \n")
        exit_program(-1)
    StdI.num_pivot = 4
    if StdI.W != 3:
        print("DEBUG: W != 3")
        exit_program(-1)
    StdI.W = 9
    fp.write("# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n")
    fp.write(f"{StdI.W} {StdI.L} {StdI.num_pivot} {StdI.ishift_nspin}\n")

    # 6-spin star list
    StdI.list_6spin_star = np.zeros((StdI.num_pivot, 7), dtype=int)

    StdI.list_6spin_star[0, :] = [1, 1, 1, 1, 4, 2, -1]
    StdI.list_6spin_star[1, :] = [6, 1, 1, 1, 6, 7, 1]
    StdI.list_6spin_star[2, :] = [6, 1, 1, 1, 4, 2, 1]
    StdI.list_6spin_star[3, :] = [5, 1, 1, 1, 4, 2, 1]

    fp.write("# StdI->list_6spin_star\n")
    for ipivot in range(StdI.num_pivot):
        fp.write(f"# pivot {ipivot}\n")
        for isite in range(7):
            fp.write(f"{StdI.list_6spin_star[ipivot, isite]} ")
        fp.write("\n")

    # 6-spin pair list
    max_kintr = max(StdI.list_6spin_star[ip, 0] for ip in range(StdI.num_pivot))
    StdI.list_6spin_pair = np.zeros((StdI.num_pivot, 7, max_kintr), dtype=int)

    # pivot 0
    StdI.list_6spin_pair[0, :, 0] = [0, 4, 1, 2, 3, 5, 3]

    # pivot 1
    StdI.list_6spin_pair[1, :, 0] = [0, 1, 2, 3, 4, 5, 3]
    StdI.list_6spin_pair[1, :, 1] = [1, 2, 0, 3, 4, 5, 1]
    StdI.list_6spin_pair[1, :, 2] = [0, 2, 1, 3, 4, 5, 2]
    StdI.list_6spin_pair[1, :, 3] = [1, 3, 0, 2, 4, 5, 3]
    StdI.list_6spin_pair[1, :, 4] = [2, 4, 0, 1, 3, 5, 2]
    StdI.list_6spin_pair[1, :, 5] = [2, 5, 0, 1, 3, 4, 1]

    # pivot 2
    StdI.list_6spin_pair[2, :, 0] = [0, 1, 2, 3, 4, 5, 3]
    StdI.list_6spin_pair[2, :, 1] = [1, 2, 0, 3, 4, 5, 1]
    StdI.list_6spin_pair[2, :, 2] = [0, 2, 1, 3, 4, 5, 2]
    StdI.list_6spin_pair[2, :, 3] = [1, 3, 0, 2, 4, 5, 3]
    StdI.list_6spin_pair[2, :, 4] = [2, 5, 0, 1, 3, 4, 2]
    StdI.list_6spin_pair[2, :, 5] = [2, 4, 0, 1, 3, 5, 1]

    # pivot 3
    StdI.list_6spin_pair[3, :, 0] = [0, 1, 2, 3, 4, 5, 3]
    StdI.list_6spin_pair[3, :, 1] = [1, 2, 0, 3, 4, 5, 1]
    StdI.list_6spin_pair[3, :, 2] = [0, 2, 1, 3, 4, 5, 2]
    StdI.list_6spin_pair[3, :, 3] = [2, 5, 0, 1, 3, 4, 2]
    StdI.list_6spin_pair[3, :, 4] = [2, 4, 0, 1, 3, 5, 1]

    fp.write("# StdI->list_6spin_pair\n")
    for ipivot in range(StdI.num_pivot):
        fp.write(f"# pivot {ipivot}\n")
        for kintr in range(StdI.list_6spin_star[ipivot, 0]):
            for isite in range(7):
                fp.write(f"{StdI.list_6spin_pair[ipivot, isite, kintr]} ")
            fp.write("\n")

    fp.close()
