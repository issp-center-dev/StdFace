"""
Standard mode for the tetragonal (square) lattice.

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


def tetragonal(StdI: StdIntList) -> None:
    """Setup a Hamiltonian for the square (tetragonal) lattice.

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

    StdI.NsiteUC = 1

    print("  @ Lattice Size & Shape\n")

    StdI.a = print_val_d("a", StdI.a, 1.0)
    StdI.length[0] = print_val_d("Wlength", StdI.length[0], StdI.a)
    StdI.length[1] = print_val_d("Llength", StdI.length[1], StdI.a)
    StdI.direct[0, 0] = print_val_d("Wx", StdI.direct[0, 0], StdI.length[0])
    StdI.direct[0, 1] = print_val_d("Wy", StdI.direct[0, 1], 0.0)
    StdI.direct[1, 0] = print_val_d("Lx", StdI.direct[1, 0], 0.0)
    StdI.direct[1, 1] = print_val_d("Ly", StdI.direct[1, 1], StdI.length[1])

    StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
    StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)

    init_site(StdI, fp, 2)
    StdI.tau[0, 0] = 0.0
    StdI.tau[0, 1] = 0.0
    StdI.tau[0, 2] = 0.0

    # (2) check & store parameters of Hamiltonian
    print("\n  @ Hamiltonian \n")
    not_used_j("J2", StdI.J2All, StdI.J2)
    not_used_j("J2'", StdI.J2pAll, StdI.J2p)
    not_used_c("t2", StdI.t2)
    not_used_d("t2'", StdI.t2p)
    not_used_d("V2", StdI.V2)
    not_used_d("V2'", StdI.V2p)
    not_used_d("K", StdI.K)
    StdI.h = print_val_d("h", StdI.h, 0.0)
    StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
    StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)

    if StdI.model == "spin":
        StdI.S2 = print_val_i("2S", StdI.S2, 1)
        StdI.D[2, 2] = print_val_d("D", StdI.D[2, 2], 0.0)
        input_spin_nn(StdI.J, StdI.JAll, StdI.J0, StdI.J0All, "J0")
        input_spin_nn(StdI.J, StdI.JAll, StdI.J1, StdI.J1All, "J1")
        input_spin_nn(StdI.Jp, StdI.JpAll, StdI.J0p, StdI.J0pAll, "J0'")
        input_spin_nn(StdI.Jp, StdI.JpAll, StdI.J1p, StdI.J1pAll, "J1'")
        input_spin_nn(StdI.Jpp, StdI.JppAll, StdI.J0pp, StdI.J0ppAll, "J0''")
        input_spin_nn(StdI.Jpp, StdI.JppAll, StdI.J1pp, StdI.J1ppAll, "J1''")

        not_used_d("mu", StdI.mu)
        not_used_d("U", StdI.U)
        not_used_c("t", StdI.t)
        not_used_c("t0", StdI.t0)
        not_used_c("t1", StdI.t1)
        not_used_c("t'", StdI.tp)
        not_used_c("t0'", StdI.t0p)
        not_used_c("t1'", StdI.t1p)
        not_used_c("t''", StdI.tpp)
        not_used_c("t0''", StdI.t0pp)
        not_used_c("t1''", StdI.t1pp)
        not_used_d("V", StdI.V)
        not_used_d("V0", StdI.V0)
        not_used_d("V1", StdI.V1)
        not_used_d("V'", StdI.Vp)
        not_used_d("V0'", StdI.V0p)
        not_used_d("V1'", StdI.V1p)
        not_used_d("V''", StdI.Vpp)
        not_used_d("V0''", StdI.V0pp)
        not_used_d("V1''", StdI.V1pp)
    else:
        StdI.mu = print_val_d("mu", StdI.mu, 0.0)
        StdI.U = print_val_d("U", StdI.U, 0.0)
        StdI.t0 = input_hopp(StdI.t, StdI.t0, "t0")
        StdI.t1 = input_hopp(StdI.t, StdI.t1, "t1")
        StdI.t0p = input_hopp(StdI.tp, StdI.t0p, "t0'")
        StdI.t1p = input_hopp(StdI.tp, StdI.t1p, "t1'")
        StdI.t0pp = input_hopp(StdI.tpp, StdI.t0pp, "t0''")
        StdI.t1pp = input_hopp(StdI.tpp, StdI.t1pp, "t1''")
        StdI.V0 = input_coulomb_v(StdI.V, StdI.V0, "V0")
        StdI.V1 = input_coulomb_v(StdI.V, StdI.V1, "V1")
        StdI.V0p = input_coulomb_v(StdI.Vp, StdI.V0p, "V0'")
        StdI.V1p = input_coulomb_v(StdI.Vp, StdI.V1p, "V1'")
        StdI.V0pp = input_coulomb_v(StdI.Vpp, StdI.V0pp, "V0''")
        StdI.V1pp = input_coulomb_v(StdI.Vpp, StdI.V1pp, "V1''")
        StdI.Vp = print_val_d("V'", StdI.Vp, 0.0)

        not_used_j("J0", StdI.J0All, StdI.J0)
        not_used_j("J1", StdI.J1All, StdI.J1)
        not_used_j("J'", StdI.JpAll, StdI.Jp)
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

    # (4) Compute upper limit of Transfer & Interaction and malloc them
    if StdI.model == "spin":
        ntransMax = StdI.nsite * (StdI.S2 + 1 + 2 * StdI.S2)
        nintrMax = (StdI.NCell * (StdI.NsiteUC + 2 + 2 + 2)
                    * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1))
    else:
        ntransMax = StdI.NCell * 2 * (2 * StdI.NsiteUC + 4 + 4 + 4)
        nintrMax = StdI.NCell * (StdI.NsiteUC + 4 * (2 + 2 + 2))
        if StdI.model == "kondo":
            ntransMax += StdI.nsite // 2 * (StdI.S2 + 1 + 2 * StdI.S2)
            nintrMax += StdI.nsite // 2 * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1)

    malloc_interactions(StdI, ntransMax, nintrMax)

    # (5) Set Transfer & Interaction
    for kCell in range(StdI.NCell):
        iW = StdI.Cell[kCell, 0]
        iL = StdI.Cell[kCell, 1]

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
                mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, jsite)

        # Nearest neighbor along W
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 1, 0, 0, 0, 1)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0, isite, jsite, dR)
            coulomb(StdI, StdI.V0, isite, jsite)

        # Nearest neighbor along L
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 0, 1, 0, 0, 1)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1, isite, jsite, dR)
            coulomb(StdI, StdI.V1, isite, jsite)

        # Second nearest neighbor W+L
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 1, 1, 0, 0, 2)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0p, isite, jsite, dR)
            coulomb(StdI, StdI.V0p, isite, jsite)

        # Second nearest neighbor W-L
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 1, -1, 0, 0, 2)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1p, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1p, isite, jsite, dR)
            coulomb(StdI, StdI.V1p, isite, jsite)

        # Third nearest neighbor along 2W
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 2, 0, 0, 0, 3)
        if StdI.model == "spin":
            general_j(StdI, StdI.J0pp, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t0pp, isite, jsite, dR)
            coulomb(StdI, StdI.V0pp, isite, jsite)

        # Third nearest neighbor along 2L
        isite, jsite, Cphase, dR = set_label(StdI, fp, iW, iL, 0, 2, 0, 0, 3)
        if StdI.model == "spin":
            general_j(StdI, StdI.J1pp, StdI.S2, StdI.S2, isite, jsite)
        else:
            hopping(StdI, Cphase * StdI.t1pp, isite, jsite, dR)
            coulomb(StdI, StdI.V1pp, isite, jsite)

    if StdI.solver != "HWAVE" or StdI.lattice_gp == 1:
        fp.write("plot '-' w d lc 7\n0.0 0.0\nend\npause -1\n")
        fp.close()

    print_geometry(StdI)
