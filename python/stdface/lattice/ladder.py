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

import logging
import numpy as np

from ..core.stdface_vals import StdIntList, ModelType
from ..core.param_check import (
    print_val_d, print_val_i,
    not_used_j, not_used_d, not_used_i, required_val_i,
)
from .input_params import input_spin, input_hopp, input_coulomb_v
from .interaction_builder import (
    malloc_interactions,
    add_neighbor_interaction, add_local_terms,
)
from .site_util import (
    init_site, set_label, set_local_spin_flags,
    new_gnuplot_buffer, GnuplotData, print_geometry,
)
from .boost_output import (
    write_boost_mag_field, write_boost_j_symmetric,
    write_boost_6spin_star, write_boost_6spin_pair,
)


logger = logging.getLogger(__name__)


def ladder(StdI: StdIntList) -> "GnuplotData | None":
    """Setup a Hamiltonian for the ladder lattice.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    buf = new_gnuplot_buffer(StdI)

    # 1. Set lattice size and shape parameters
    logger.info("  @ Lattice Size & Shape\n")

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
    init_site(StdI, 2)

    for isite in range(StdI.NsiteUC):
        StdI.tau[isite, 0] = float(isite) / float(StdI.NsiteUC)
        StdI.tau[isite, 1] = 0.0
        StdI.tau[isite, 2] = 0.0

    # 2. Set Hamiltonian parameters
    logger.info("\n  @ Hamiltonian \n")

    not_used_j("J", StdI.JAll, StdI.J)
    not_used_j("J'", StdI.JpAll, StdI.Jp)
    not_used_d("t", StdI.t)
    not_used_d("t'", StdI.tp)
    not_used_d("V", StdI.V)
    not_used_d("V'", StdI.Vp)
    not_used_d("K", StdI.K)

    StdI.h = print_val_d("h", StdI.h, 0.0)
    StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
    StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)

    if StdI.model == ModelType.SPIN:
        StdI.S2 = print_val_i("2S", StdI.S2, 1)
        StdI.D[2, 2] = print_val_d("D", StdI.D[2, 2], 0.0)
        input_spin(StdI.J0, StdI.J0All, "J0")
        input_spin(StdI.J1, StdI.J1All, "J1")
        input_spin(StdI.J2, StdI.J2All, "J2")
        input_spin(StdI.J1p, StdI.J1pAll, "J1'")
        input_spin(StdI.J2p, StdI.J2pAll, "J2'")

        not_used_d("mu", StdI.mu)
        not_used_d("U", StdI.U)
        not_used_d("t0", StdI.t0)
        not_used_d("t1", StdI.t1)
        not_used_d("t2", StdI.t2)
        not_used_d("t1'", StdI.t1p)
        not_used_d("t2'", StdI.t2p)
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

        if StdI.model == ModelType.HUBBARD:
            not_used_i("2S", StdI.S2)
            not_used_j("J", StdI.JAll, StdI.J)
        else:
            StdI.S2 = print_val_i("2S", StdI.S2, 1)
            input_spin(StdI.J, StdI.JAll, "J")

    logger.info("\n  @ Numerical conditions\n")

    # 3. Set local spin flags and number of sites
    set_local_spin_flags(StdI, StdI.L * StdI.NsiteUC)

    # 4. Calculate maximum number of transfer terms (for pump arrays) and allocate
    if StdI.model == ModelType.SPIN:
        ntransMax = StdI.L * StdI.NsiteUC * (StdI.S2 + 1 + 2 * StdI.S2)
    else:
        ntransMax = (StdI.L * StdI.NsiteUC * 2 * (2 + 2 + 2)
                     + StdI.L * (StdI.NsiteUC - 1) * 2 * (2 + 2 + 2))
        if StdI.model == ModelType.KONDO:
            ntransMax += StdI.L * StdI.NsiteUC * (StdI.S2 + 1 + 2 * StdI.S2)

    malloc_interactions(StdI, ntransMax)

    # 5. Set all interactions
    for cell_l in range(StdI.L):
        for uc_i in range(StdI.NsiteUC):
            isite = uc_i + cell_l * StdI.NsiteUC
            if StdI.model == ModelType.KONDO:
                isite += StdI.L * StdI.NsiteUC

            # Local terms
            add_local_terms(StdI, isite, uc_i + cell_l * StdI.NsiteUC)

            # Leg bonds: (dW, dL, sj_offset, nn, J, t, V)
            _LEG_BONDS = (
                (0, 1, 0, 1, StdI.J1, StdI.t1, StdI.V1),   # nn along ladder
                (0, 2, 0, 2, StdI.J1p, StdI.t1p, StdI.V1p), # nnn along ladder
            )
            for dW, dL, sj_off, nn, J, t, V in _LEG_BONDS:
                add_neighbor_interaction(
                    StdI, buf, 0, cell_l, dW, dL, uc_i, uc_i + sj_off, nn, J, t, V)

            # Rung/diagonal bonds (only between adjacent legs)
            if uc_i < StdI.NsiteUC - 1:
                _RUNG_BONDS = (
                    (0, 0, 1, StdI.J0, StdI.t0, StdI.V0),    # vertical
                    (0, 1, 1, StdI.J2, StdI.t2, StdI.V2),    # diagonal 1
                    (0, -1, 1, StdI.J2p, StdI.t2p, StdI.V2p), # diagonal 2
                )
                for dW, dL, nn, J, t, V in _RUNG_BONDS:
                    add_neighbor_interaction(
                        StdI, buf, 0, cell_l, dW, dL, uc_i, uc_i + 1, nn, J, t, V)
    return buf.build(StdI) if buf else None


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

    with open("boost.def", "w") as fp:
        # Magnetic field
        write_boost_mag_field(fp, StdI)

        # Interaction parameters
        fp.write(f"{5}  # Number of type of J\n")

        # J1 - Vertical interactions
        fp.write("# J 1 (inter chain, vertical)\n")
        write_boost_j_symmetric(fp, StdI.J0)

        # J2 - Nearest neighbor along chain
        fp.write("# J 2 (Nearest neighbor, along chain)\n")
        write_boost_j_symmetric(fp, StdI.J1)

        # J3 - Second nearest neighbor along chain
        fp.write("# J 3 (Second nearest neighbor, along chain)\n")
        write_boost_j_symmetric(fp, StdI.J1p)

        # J4 - Diagonal 1
        fp.write("# J 4 (inter chain, diagonal1)\n")
        write_boost_j_symmetric(fp, StdI.J2)

        # J5 - Diagonal 2
        fp.write("# J 5 (inter chain, diagonal2)\n")
        write_boost_j_symmetric(fp, StdI.J2p)

        # Validate parameters
        if StdI.S2 != 1:
            msg = "\n ERROR! S2 must be 1 in Boost. \n"
            logger.error(msg)
            raise ValueError(msg)
        StdI.ishift_nspin = 2
        if StdI.W != 2:
            msg = "\n ERROR! W != 2 \n"
            logger.error(msg)
            raise ValueError(msg)
        if StdI.L % 2 != 0:
            msg = "\n ERROR! L %% 2 != 0 \n"
            logger.error(msg)
            raise ValueError(msg)
        if StdI.L < 4:
            msg = "\n ERROR! L < 4 \n"
            logger.error(msg)
            raise ValueError(msg)

        StdI.W = StdI.L
        StdI.L = 2
        StdI.num_pivot = StdI.W // 2

        fp.write("# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n")
        fp.write(f"{StdI.W} {StdI.L} {StdI.num_pivot} {StdI.ishift_nspin}\n")

        # 6-spin star list
        StdI.list_6spin_star = np.zeros((StdI.num_pivot, 7), dtype=int)
        for ipivot in range(StdI.num_pivot):
            StdI.list_6spin_star[ipivot, :] = [7, 1, 1, 1, 1, 1, 1]

        write_boost_6spin_star(fp, StdI)

        # 6-spin pair list
        StdI.list_6spin_pair = np.zeros((StdI.num_pivot, 7, 7), dtype=int)
        for ipivot in range(StdI.num_pivot):
            StdI.list_6spin_pair[ipivot, :, 0] = [0, 1, 2, 3, 4, 5, 1]
            StdI.list_6spin_pair[ipivot, :, 1] = [0, 2, 1, 3, 4, 5, 2]
            StdI.list_6spin_pair[ipivot, :, 2] = [1, 3, 0, 2, 4, 5, 2]
            StdI.list_6spin_pair[ipivot, :, 3] = [0, 4, 1, 2, 3, 5, 3]
            StdI.list_6spin_pair[ipivot, :, 4] = [1, 5, 0, 2, 3, 4, 3]
            StdI.list_6spin_pair[ipivot, :, 5] = [0, 3, 1, 2, 4, 5, 4]
            StdI.list_6spin_pair[ipivot, :, 6] = [1, 2, 0, 3, 4, 5, 5]

        write_boost_6spin_pair(fp, StdI)


# ---------------------------------------------------------------------------
#  Lattice plugin registration
# ---------------------------------------------------------------------------

from . import LatticePlugin, register_lattice

_ladder_setup = ladder
_ladder_boost = ladder_boost


class LadderPlugin(LatticePlugin):
    """Plugin for the 1D ladder lattice."""

    @property
    def name(self) -> str:
        return "ladder"

    @property
    def aliases(self) -> list[str]:
        return ["ladder", "ladderlattice"]

    @property
    def ndim(self) -> int:
        return 1

    def setup(self, StdI: StdIntList) -> "GnuplotData | None":
        """Delegate to the ladder() function."""
        return _ladder_setup(StdI)

    def boost(self, StdI: StdIntList) -> None:
        """Delegate to the ladder_boost() function."""
        _ladder_boost(StdI)


register_lattice(LadderPlugin())
