"""
Standard mode for the face-centered orthorhombic lattice.

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
    print_val_d,
    print_val_i,
    init_site,
    find_site,
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
    print_xsf,
)


def fc_ortho(StdI: StdIntList) -> None:
    """Setup a Hamiltonian for the face-centered orthorhombic lattice.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    # (1) Compute the shape of the super-cell and sites in the super-cell
    fp = open("lattice.xsf", "w")

    StdI.NsiteUC = 1

    print("  @ Lattice Size & Shape\n")

    StdI.a = print_val_d("a", StdI.a, 1.0)
    StdI.length[0] = print_val_d("Wlength", StdI.length[0], StdI.a)
    StdI.length[1] = print_val_d("Llength", StdI.length[1], StdI.a)
    StdI.length[2] = print_val_d("Hlength", StdI.length[2], StdI.a)
    StdI.direct[0, 0] = print_val_d("Wx", StdI.direct[0, 0], 0.0)
    StdI.direct[0, 1] = print_val_d("Wy", StdI.direct[0, 1], 0.5 * StdI.length[1])
    StdI.direct[0, 2] = print_val_d("Wz", StdI.direct[0, 2], 0.5 * StdI.length[2])
    StdI.direct[1, 0] = print_val_d("Lx", StdI.direct[1, 0], 0.5 * StdI.length[0])
    StdI.direct[1, 1] = print_val_d("Ly", StdI.direct[1, 1], 0.0)
    StdI.direct[1, 2] = print_val_d("Lz", StdI.direct[1, 2], 0.5 * StdI.length[2])
    StdI.direct[2, 0] = print_val_d("Hx", StdI.direct[2, 0], 0.5 * StdI.length[0])
    StdI.direct[2, 1] = print_val_d("Hy", StdI.direct[2, 1], 0.5 * StdI.length[1])
    StdI.direct[2, 2] = print_val_d("Hz", StdI.direct[2, 2], 0.0)

    StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
    StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)
    StdI.phase[2] = print_val_d("phase2", StdI.phase[2], 0.0)

    init_site(StdI, fp, 3)
    StdI.tau[0, 0] = 0.0
    StdI.tau[0, 1] = 0.0
    StdI.tau[0, 2] = 0.0

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
        input_spin_nn(StdI.Jpp, StdI.JppAll, StdI.J0pp, StdI.J0ppAll, "J0''")
        input_spin_nn(StdI.Jpp, StdI.JppAll, StdI.J1pp, StdI.J1ppAll, "J1''")
        input_spin_nn(StdI.Jpp, StdI.JppAll, StdI.J2pp, StdI.J2ppAll, "J2''")

        not_used_d("mu", StdI.mu)
        not_used_d("U", StdI.U)
        not_used_c("t", StdI.t)
        not_used_c("t0", StdI.t0)
        not_used_c("t1", StdI.t1)
        not_used_c("t2", StdI.t2)
        not_used_c("t'", StdI.tp)
        not_used_c("t0'", StdI.t0p)
        not_used_c("t1'", StdI.t1p)
        not_used_c("t2'", StdI.t2p)
        not_used_c("t''", StdI.tpp)
        not_used_d("V", StdI.V)
        not_used_d("V0", StdI.V0)
        not_used_d("V1", StdI.V1)
        not_used_d("V'", StdI.Vp)
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
        not_used_j("J0'", StdI.J0pAll, StdI.J0p)
        not_used_j("J1'", StdI.J1pAll, StdI.J1p)
        not_used_j("J2'", StdI.J2pAll, StdI.J2p)
        not_used_j("J''", StdI.JppAll, StdI.Jpp)
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
        nintrMax = (StdI.NCell * (StdI.NsiteUC + 6 + 3 + 0)
                    * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1))
    else:
        ntransMax = StdI.NCell * 2 * (2 * StdI.NsiteUC + 12 + 6 + 0)
        nintrMax = StdI.NCell * (StdI.NsiteUC + 4 * (6 + 3 + 0))
        if StdI.model == "kondo":
            ntransMax += StdI.nsite // 2 * (StdI.S2 + 1 + 2 * StdI.S2)
            nintrMax += StdI.nsite // 2 * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1)

    malloc_interactions(StdI, ntransMax, nintrMax)

    # (5) Set Transfer & Interaction
    for kCell in range(StdI.NCell):
        iW = StdI.Cell[kCell, 0]
        iL = StdI.Cell[kCell, 1]
        iH = StdI.Cell[kCell, 2]

        # Local term
        isite = kCell
        if StdI.model == "kondo":
            isite += StdI.NCell

        if StdI.model == "spin":
            mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, isite)
            general_j(StdI, StdI.D, StdI.S2, StdI.S2, isite, isite)
        else:
            hubbard_local(StdI, StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, StdI.U, isite)
            if StdI.model == "kondo":
                jsite = kCell
                general_j(StdI, StdI.J, 1, StdI.S2, isite, jsite)

        # (2) Nearest neighbor along W
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, 1, 0, 0, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0, isite, jsite, dR)
            coulomb(StdI, StdI.V0, isite, jsite)

        # (3) Nearest neighbor along W (equivalent direction)
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, 0, 1, -1, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0, isite, jsite, dR)
            coulomb(StdI, StdI.V0, isite, jsite)

        # (4) Nearest neighbor along L
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, 0, 1, 0, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1, isite, jsite, dR)
            coulomb(StdI, StdI.V1, isite, jsite)

        # (5) Nearest neighbor along L (equivalent direction)
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, -1, 0, 1, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1, isite, jsite, dR)
            coulomb(StdI, StdI.V1, isite, jsite)

        # (6) Nearest neighbor along H
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, 0, 0, 1, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J2, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t2, isite, jsite, dR)
            coulomb(StdI, StdI.V2, isite, jsite)

        # (7) Nearest neighbor along H (equivalent direction)
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, 1, -1, 0, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J2, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t2, isite, jsite, dR)
            coulomb(StdI, StdI.V2, isite, jsite)

        # (8) Second nearest neighbor along -W+L+H
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, -1, 1, 1, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0p, isite, jsite, dR)
            coulomb(StdI, StdI.V0p, isite, jsite)

        # (9) Second nearest neighbor along -L+H+W
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, 1, -1, 1, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1p, isite, jsite, dR)
            coulomb(StdI, StdI.V1p, isite, jsite)

        # (10) Second nearest neighbor along -H+W+L
        isite, jsite, Cphase, dR = find_site(StdI, iW, iL, iH, 1, 1, -1, 0, 0)
        if StdI.model == "spin":
            general_j(StdI, StdI.J2p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t2p, isite, jsite, dR)
            coulomb(StdI, StdI.V2p, isite, jsite)

    fp.close()
    print_xsf(StdI)
    print_geometry(StdI)
