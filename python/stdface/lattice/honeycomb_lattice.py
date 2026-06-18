"""
Standard mode for the honeycomb lattice.

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

import logging
import math

import numpy as np

from ..core.stdface_vals import StdIntList, ModelType
from ..core.param_check import (
    print_val_d, print_val_i,
    not_used_j, not_used_d, not_used_i,
)
from .input_params import input_spin_nn, input_spin, input_hopp, input_coulomb_v
from .interaction_builder import (
    compute_max_interactions, malloc_interactions,
    add_neighbor_interaction, add_local_terms,
)
from .site_util import (
    init_site, set_label, set_local_spin_flags,
    lattice_gp,
)
from .boost_output import (
    write_boost_mag_field, write_boost_j_symmetric,
    write_boost_6spin_star, write_boost_6spin_pair,
)


logger = logging.getLogger(__name__)


def honeycomb(StdI: StdIntList) -> None:
    """Setup a Hamiltonian for the honeycomb lattice.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    # (1) Compute the shape of the super-cell and sites in the super-cell
    with lattice_gp(StdI) as fp:

        StdI.NsiteUC = 2

        logger.info("  @ Lattice Size & Shape\n")

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
        StdI.tau[1, 0] = 1.0 / 3.0
        StdI.tau[1, 1] = 1.0 / 3.0
        StdI.tau[1, 2] = 0.0

        # (2) check & store parameters of Hamiltonian
        logger.info("\n  @ Hamiltonian \n")
        not_used_d("K", StdI.K)
        StdI.h = print_val_d("h", StdI.h, 0.0)
        StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
        StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)

        if StdI.model == ModelType.SPIN:
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
            not_used_d("t", StdI.t)
            not_used_d("t0", StdI.t0)
            not_used_d("t1", StdI.t1)
            not_used_d("t2", StdI.t2)
            not_used_d("t'", StdI.tp)
            not_used_d("t0'", StdI.t0p)
            not_used_d("t1'", StdI.t1p)
            not_used_d("t2'", StdI.t2p)
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
            StdI.t0pp = input_hopp(StdI.tpp, StdI.t0pp, "t0''")
            StdI.t1pp = input_hopp(StdI.tpp, StdI.t1pp, "t1''")
            StdI.t2pp = input_hopp(StdI.tpp, StdI.t2pp, "t2''")
            StdI.V0 = input_coulomb_v(StdI.V, StdI.V0, "V0")
            StdI.V1 = input_coulomb_v(StdI.V, StdI.V1, "V1")
            StdI.V2 = input_coulomb_v(StdI.V, StdI.V2, "V2")
            StdI.V0p = input_coulomb_v(StdI.Vp, StdI.V0p, "V0'")
            StdI.V1p = input_coulomb_v(StdI.Vp, StdI.V1p, "V1'")
            StdI.V2p = input_coulomb_v(StdI.Vp, StdI.V2p, "V2'")
            StdI.V0pp = input_coulomb_v(StdI.Vpp, StdI.V0pp, "V0''")
            StdI.V1pp = input_coulomb_v(StdI.Vpp, StdI.V1pp, "V1''")
            StdI.V2pp = input_coulomb_v(StdI.Vpp, StdI.V2pp, "V2''")
            StdI.Vp = print_val_d("V'", StdI.Vp, 0.0)

            not_used_j("J0", StdI.J0All, StdI.J0)
            not_used_j("J1", StdI.J1All, StdI.J1)
            not_used_j("J2", StdI.J2All, StdI.J2)
            not_used_j("J'", StdI.JpAll, StdI.Jp)
            not_used_d("D", StdI.D[2, 2])

            if StdI.model == ModelType.HUBBARD:
                not_used_i("2S", StdI.S2)
                not_used_j("J", StdI.JAll, StdI.J)
            else:
                StdI.S2 = print_val_i("2S", StdI.S2, 1)
                input_spin(StdI.J, StdI.JAll, "J")

        logger.info("\n  @ Numerical conditions\n")

        # (3) Set local spin flag and number of sites
        set_local_spin_flags(StdI, StdI.NsiteUC * StdI.NCell)

        # (4) Compute upper limit of Transfer & Interaction
        #     nn=3, nnn=6, nnnn=3 → n_bonds=12
        ntransMax, nintrMax = compute_max_interactions(StdI, n_bonds=3 + 6 + 3)
        malloc_interactions(StdI, ntransMax, nintrMax)

        # (5) Set Transfer & Interaction
        for kCell in range(StdI.NCell):
            cell_w = StdI.Cell[kCell, 0]
            cell_l = StdI.Cell[kCell, 1]

            # Local term
            isite = StdI.NsiteUC * kCell
            if StdI.model == ModelType.KONDO:
                isite += StdI.NsiteUC * StdI.NCell
            jsite_base = StdI.NsiteUC * kCell
            for uc_i in range(StdI.NsiteUC):
                add_local_terms(StdI, isite + uc_i, jsite_base + uc_i)

            # Neighbor bonds: (dW, dL, site_i, site_j, nn_level, J, t, V)
            _BONDS = (
                # Nearest neighbor (nn=1)
                (0, 0, 0, 1, 1, StdI.J0, StdI.t0, StdI.V0),     # intra cell
                (1, 0, 1, 0, 1, StdI.J1, StdI.t1, StdI.V1),     # along W
                (0, 1, 1, 0, 1, StdI.J2, StdI.t2, StdI.V2),     # along L
                # Second nearest neighbor (nn=2)
                (1, 0, 0, 0, 2, StdI.J2p, StdI.t2p, StdI.V2p),  # along W, 0->0
                (1, 0, 1, 1, 2, StdI.J2p, StdI.t2p, StdI.V2p),  # along W, 1->1
                (0, 1, 0, 0, 2, StdI.J1p, StdI.t1p, StdI.V1p),  # along L, 0->0
                (0, 1, 1, 1, 2, StdI.J1p, StdI.t1p, StdI.V1p),  # along L, 1->1
                (1, -1, 0, 0, 2, StdI.J0p, StdI.t0p, StdI.V0p), # along W-L, 0->0
                (1, -1, 1, 1, 2, StdI.J0p, StdI.t0p, StdI.V0p), # along W-L, 1->1
                # Third nearest neighbor (nn=3)
                (1, -1, 0, 1, 3, StdI.J1pp, StdI.t1pp, StdI.V1pp),  # along W-L
                (-1, -1, 0, 1, 3, StdI.J0pp, StdI.t0pp, StdI.V0pp), # along -W-L
                (-1, 1, 0, 1, 3, StdI.J2pp, StdI.t2pp, StdI.V2pp),  # along -W+L
            )
            for dW, dL, si, sj, nn, J, t, V in _BONDS:
                add_neighbor_interaction(
                    StdI, fp, cell_w, cell_l, dW, dL, si, sj, nn, J, t, V)


def honeycomb_boost(StdI: StdIntList) -> None:
    """Setup a boosted Hamiltonian for the honeycomb lattice (HPhi only).

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    if StdI.box[0, 1] != 0 or StdI.box[1, 0] != 0:
        msg = "\nERROR ! (a0W, a0L, a1W, a1L) can not be used with SpinGCBoost.\n"
        logger.error(msg)
        raise ValueError(msg)

    if np.any(np.abs(StdI.Jp) > 1.0e-8):
        msg = "\nERROR ! J' can not be used with SpinGCBoost.\n"
        logger.error(msg)
        raise ValueError(msg)

    # Magnetic field
    with open("boost.def", "w") as fp:
        write_boost_mag_field(fp, StdI)

        # Interaction
        fp.write(f"{3}  # Number of type of J\n")
        fp.write("# J 0\n")
        write_boost_j_symmetric(fp, StdI.J0)
        fp.write("# J 1\n")
        write_boost_j_symmetric(fp, StdI.J1)
        fp.write("# J 2\n")
        write_boost_j_symmetric(fp, StdI.J2)

        # Topology
        if StdI.S2 != 1:
            msg = "\n ERROR! S2 must be 1 in Boost. \n"
            logger.error(msg)
            raise ValueError(msg)
        StdI.ishift_nspin = 3
        if StdI.L < 2:
            msg = "\n ERROR! L < 2 \n"
            logger.error(msg)
            raise ValueError(msg)
        if StdI.W % StdI.ishift_nspin != 0:
            msg = f"\n ERROR! W %% {StdI.ishift_nspin} != 0 \n"
            logger.error(msg)
            raise ValueError(msg)
        StdI.num_pivot = 2
        if StdI.W != 3:
            msg = "DEBUG: W != 3"
            logger.error(msg)
            raise ValueError(msg)
        StdI.W = 6
        fp.write("# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n")
        fp.write(f"{StdI.W} {StdI.L} {StdI.num_pivot} {StdI.ishift_nspin}\n")

        # 6-spin star list
        StdI.list_6spin_star = np.zeros((StdI.num_pivot, 7), dtype=int)

        StdI.list_6spin_star[0, :] = [5, 1, 1, 1, 2, 1, 1]
        StdI.list_6spin_star[1, :] = [4, 1, 1, 1, 2, 2, 1]

        write_boost_6spin_star(fp, StdI)

        # 6-spin pair list
        max_kintr = max(StdI.list_6spin_star[ip, 0] for ip in range(StdI.num_pivot))
        StdI.list_6spin_pair = np.zeros((StdI.num_pivot, 7, max_kintr), dtype=int)

        # pivot 0
        StdI.list_6spin_pair[0, :, 0] = [0, 1, 2, 3, 4, 5, 1]
        StdI.list_6spin_pair[0, :, 1] = [1, 2, 0, 3, 4, 5, 2]
        StdI.list_6spin_pair[0, :, 2] = [2, 3, 0, 1, 4, 5, 1]
        StdI.list_6spin_pair[0, :, 3] = [0, 4, 1, 2, 3, 5, 2]
        StdI.list_6spin_pair[0, :, 4] = [1, 5, 0, 2, 3, 4, 3]

        # pivot 1
        StdI.list_6spin_pair[1, :, 0] = [0, 1, 2, 3, 4, 5, 2]
        StdI.list_6spin_pair[1, :, 1] = [1, 2, 0, 3, 4, 5, 1]
        StdI.list_6spin_pair[1, :, 2] = [0, 4, 1, 2, 3, 5, 3]
        StdI.list_6spin_pair[1, :, 3] = [2, 5, 0, 1, 3, 4, 3]

        write_boost_6spin_pair(fp, StdI)


# ---------------------------------------------------------------------------
#  Lattice plugin registration
# ---------------------------------------------------------------------------

from . import LatticePlugin, register_lattice


class HoneycombPlugin(LatticePlugin):
    """Plugin for the 2D honeycomb lattice."""

    @property
    def name(self) -> str:
        return "honeycomb"

    @property
    def aliases(self) -> list[str]:
        return ["honeycomb", "honeycomblattice"]

    @property
    def ndim(self) -> int:
        return 2

    def setup(self, StdI: StdIntList) -> None:
        """Delegate to the honeycomb() function."""
        honeycomb(StdI)

    def boost(self, StdI: StdIntList) -> None:
        """Delegate to the honeycomb_boost() function."""
        honeycomb_boost(StdI)


register_lattice(HoneycombPlugin())
