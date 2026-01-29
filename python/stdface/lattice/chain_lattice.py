"""
Standard mode for the chain lattice.

This module sets up the Hamiltonian for a 1D chain lattice model,
supporting spin, Hubbard, and Kondo models.

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

from ..core.stdface_vals import StdIntList, ModelType, SolverType
from ..core.param_check import (
    exit_program, print_val_d, print_val_i,
    not_used_d, not_used_i, not_used_j, required_val_i,
)
from .input_params import input_spin_nn, input_spin, input_hopp, input_coulomb_v
from .interaction_builder import (
    malloc_interactions,
    add_neighbor_interaction, add_local_terms,
)
from .site_util import (
    init_site, set_label, set_local_spin_flags,
    lattice_gp,
)
from .boost_output import (
    write_boost_mag_field, write_boost_j_full,
    write_boost_6spin_star, write_boost_6spin_pair,
)


def chain(StdI: StdIntList) -> None:
    """Set up a Hamiltonian for the Hubbard model on a chain lattice.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.

    Notes
    -----
    This function handles three different model types:

    - Spin model
    - Hubbard model
    - Kondo model

    The function performs the following steps:

    1. Computes super-cell shape and sites
    2. Validates and stores Hamiltonian parameters
    3. Sets local spin flags and number of sites
    4. Allocates memory for interactions
    5. Sets up transfers and interactions between sites

    The lattice geometry is written to ``lattice.gp`` for visualization.
    For the Kondo model, the number of sites is doubled to account for
    localized spins.
    """
    # (1) Compute the shape of the super-cell and sites in the super-cell
    with lattice_gp(StdI) as fp:

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
        not_used_d("phase1", StdI.phase[1])
        StdI.phase[1] = StdI.phase[0]
        StdI.phase[0] = 0.0

        required_val_i("L", StdI.L)
        not_used_i("W", StdI.W)
        StdI.W = 1
        init_site(StdI, fp, 2)
        StdI.tau[0, 0] = 0.0
        StdI.tau[0, 1] = 0.0
        StdI.tau[0, 2] = 0.0

        # (2) Check & store parameters of Hamiltonian
        print("\n  @ Hamiltonian \n")
        not_used_j("J1", StdI.J1All, StdI.J1)
        not_used_j("J2", StdI.J2All, StdI.J2)
        not_used_j("J1'", StdI.J1pAll, StdI.J1p)
        not_used_j("J2'", StdI.J2pAll, StdI.J2p)
        not_used_d("t1", StdI.t1)
        not_used_d("t2", StdI.t2)
        not_used_d("t1'", StdI.t1p)
        not_used_d("t2'", StdI.t2p)
        not_used_d("V1", StdI.V1)
        not_used_d("V2", StdI.V2)
        not_used_d("V1'", StdI.V1p)
        not_used_d("V2'", StdI.V2p)
        not_used_d("K", StdI.K)
        StdI.h = print_val_d("h", StdI.h, 0.0)
        StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
        StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)

        if StdI.model == ModelType.SPIN:
            StdI.S2 = print_val_i("2S", StdI.S2, 1)
            StdI.D[2, 2] = print_val_d("D", StdI.D[2, 2], 0.0)
            input_spin_nn(StdI.J, StdI.JAll, StdI.J0, StdI.J0All, "J0")
            input_spin_nn(StdI.Jp, StdI.JpAll, StdI.J0p, StdI.J0pAll, "J0'")
            input_spin_nn(StdI.Jpp, StdI.JppAll, StdI.J0pp, StdI.J0ppAll, "J0''")

            not_used_d("mu", StdI.mu)
            not_used_d("U", StdI.U)
            not_used_d("t", StdI.t)
            not_used_d("t0", StdI.t0)
            not_used_d("t'", StdI.tp)
            not_used_d("V", StdI.V)
            not_used_d("V0", StdI.V0)
            not_used_d("V'", StdI.Vp)
        else:
            StdI.mu = print_val_d("mu", StdI.mu, 0.0)
            StdI.U = print_val_d("U", StdI.U, 0.0)
            StdI.t0 = input_hopp(StdI.t, StdI.t0, "t0")
            StdI.t0p = input_hopp(StdI.tp, StdI.t0p, "t0'")
            StdI.t0pp = input_hopp(StdI.tpp, StdI.t0pp, "t0''")
            StdI.V0 = input_coulomb_v(StdI.V, StdI.V0, "V0")
            StdI.V0p = input_coulomb_v(StdI.Vp, StdI.V0p, "V0'")
            StdI.V0pp = input_coulomb_v(StdI.Vpp, StdI.V0pp, "V0''")

            not_used_j("J0", StdI.J0All, StdI.J0)
            not_used_j("J0'", StdI.J0pAll, StdI.J0p)
            not_used_j("J0''", StdI.J0ppAll, StdI.J0pp)
            not_used_d("D", StdI.D[2, 2])

            if StdI.model == ModelType.HUBBARD:
                not_used_i("2S", StdI.S2)
                not_used_j("J", StdI.JAll, StdI.J)
            elif StdI.model == ModelType.KONDO:
                StdI.S2 = print_val_i("2S", StdI.S2, 1)
                input_spin(StdI.J, StdI.JAll, "J")

        print("\n  @ Numerical conditions\n")

        # (3) Set local spin flag and the number of sites
        set_local_spin_flags(StdI, StdI.L)

        # (4) Compute upper limit of Transfer & Interaction and allocate them
        if StdI.model == ModelType.SPIN:
            ntransMax = StdI.L * (StdI.S2 + 1 + 2 * StdI.S2)
            nintrMax = (StdI.L * (StdI.NsiteUC + 1 + 1 + 1)
                        * (3 * StdI.S2 + 1) * (3 * StdI.S2 + 1))
        else:
            ntransMax = StdI.L * 2 * (2 * StdI.NsiteUC + 2 + 2 + 2)
            nintrMax = StdI.L * (StdI.NsiteUC + 4 * (1 + 1 + 1))

            if StdI.model == ModelType.KONDO:
                ntransMax += StdI.L * (StdI.S2 + 1 + 2 * StdI.S2)
                nintrMax += StdI.nsite // 2 * (3 * 1 + 1) * (3 * StdI.S2 + 1)

        malloc_interactions(StdI, ntransMax, nintrMax)

        # (5) Set Transfer & Interaction
        for iL in range(StdI.L):

            isite = iL
            if StdI.model == ModelType.KONDO:
                isite += StdI.L

            # Local term
            add_local_terms(StdI, isite, iL)

            # Neighbor bonds: (dW, dL, site_i, site_j, nn_level, J, t, V)
            _BONDS = (
                (0, 1, 0, 0, 1, StdI.J0, StdI.t0, StdI.V0),       # nn
                (0, 2, 0, 0, 2, StdI.J0p, StdI.t0p, StdI.V0p),    # nnn
                (0, 3, 0, 0, 3, StdI.J0pp, StdI.t0pp, StdI.V0pp), # nnnn
            )
            for dW, dL, si, sj, nn, J, t, V in _BONDS:
                add_neighbor_interaction(
                    StdI, fp, 0, iL, dW, dL, si, sj, nn, J, t, V)


def chain_boost(StdI: StdIntList) -> None:
    """Set up a Hamiltonian for the generalized Heisenberg model on a chain
    lattice using HPhi boost mode.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.

    Notes
    -----
    This function handles:

    - Magnetic field terms
    - Nearest and next-nearest neighbor interactions
    - Specialized 6-spin interactions

    The function performs:

    1. Sets up unit cell and lattice parameters
    2. Writes magnetic field configuration to ``boost.def``
    3. Writes interaction parameters
    4. Sets up topology and pivot sites
    5. Configures 6-spin interaction lists

    Only available when ``StdI.solver == "HPhi"``.

    Warnings
    --------
    - ``S2`` must be 1 in Boost mode.
    - ``L`` must be divisible by 8.
    """
    if StdI.solver != SolverType.HPhi:
        return

    StdI.NsiteUC = 1

    # Magnetic field
    with open("boost.def", "w") as fp:
        write_boost_mag_field(fp, StdI)

        # Interaction
        fp.write(f"{2}  # Number of type of J\n")
        fp.write("# J 1\n")
        write_boost_j_full(fp, StdI.J0)
        fp.write("# J 2\n")
        write_boost_j_full(fp, StdI.Jp)

        # Topology
        if StdI.S2 != 1:
            print("\n ERROR! S2 must be 1 in Boost. \n")
            exit_program(-1)
        StdI.ishift_nspin = 4
        if StdI.L % 8 != 0:
            print("\n ERROR! L % 8 != 0 \n")
            exit_program(-1)
        StdI.W = StdI.L // 2
        StdI.L = 2
        StdI.num_pivot = StdI.W // 4

        fp.write("# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n")
        fp.write(f"{StdI.W} {StdI.L} {StdI.num_pivot} {StdI.ishift_nspin}\n")

        # list_6spin_star: shape (num_pivot, 7)
        StdI.list_6spin_star = np.zeros((StdI.num_pivot, 7), dtype=int)
        for ipivot in range(StdI.num_pivot):
            StdI.list_6spin_star[ipivot, :] = [8, 1, 1, 1, 1, 1, 1]

        write_boost_6spin_star(fp, StdI)

        # list_6spin_pair: shape (num_pivot, 7, 8)
        StdI.list_6spin_pair = np.zeros((StdI.num_pivot, 7, 8), dtype=int)
        for ipivot in range(StdI.num_pivot):
            StdI.list_6spin_pair[ipivot, :, 0] = [0, 1, 2, 3, 4, 5, 1]
            StdI.list_6spin_pair[ipivot, :, 1] = [1, 2, 0, 3, 4, 5, 1]
            StdI.list_6spin_pair[ipivot, :, 2] = [2, 3, 0, 1, 4, 5, 1]
            StdI.list_6spin_pair[ipivot, :, 3] = [3, 4, 0, 1, 2, 5, 1]
            StdI.list_6spin_pair[ipivot, :, 4] = [0, 2, 1, 3, 4, 5, 2]
            StdI.list_6spin_pair[ipivot, :, 5] = [1, 3, 0, 2, 4, 5, 2]
            StdI.list_6spin_pair[ipivot, :, 6] = [2, 4, 0, 1, 3, 5, 2]
            StdI.list_6spin_pair[ipivot, :, 7] = [3, 5, 0, 1, 2, 4, 2]

        write_boost_6spin_pair(fp, StdI)


# ---------------------------------------------------------------------------
#  Lattice plugin registration
# ---------------------------------------------------------------------------

from . import LatticePlugin, register_lattice


class ChainPlugin(LatticePlugin):
    """Plugin for the 1D chain lattice."""

    @property
    def name(self) -> str:
        return "chain"

    @property
    def aliases(self) -> list[str]:
        return ["chain", "chainlattice"]

    @property
    def ndim(self) -> int:
        return 1

    def setup(self, StdI: StdIntList) -> None:
        """Delegate to the chain() function."""
        chain(StdI)

    def boost(self, StdI: StdIntList) -> None:
        """Delegate to the chain_boost() function."""
        chain_boost(StdI)


register_lattice(ChainPlugin())
