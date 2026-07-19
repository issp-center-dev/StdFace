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

import logging

from ..core.stdface_vals import StdIntList, ModelType
from ..core.param_check import (
    print_val_d, print_val_i,
    not_used_j, not_used_d, not_used_i,
)
from .input_params import input_spin_nn, input_spin, input_hopp, input_coulomb_v
from .interaction_builder import (
    compute_max_interactions, malloc_interactions,
    expand_bonds_2d,
)
from .site_util import (
    init_site, set_local_spin_flags,
    new_gnuplot_buffer, GnuplotData,
)


logger = logging.getLogger(__name__)


def tetragonal(StdI: StdIntList) -> "GnuplotData | None":
    """Setup a Hamiltonian for the square (tetragonal) lattice.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    # (1) Compute the shape of the super-cell and sites in the super-cell
    buf = new_gnuplot_buffer(StdI)

    StdI.NsiteUC = 1

    logger.info("  @ Lattice Size & Shape\n")

    StdI.a = print_val_d("a", StdI.a, 1.0)
    StdI.length[0] = print_val_d("Wlength", StdI.length[0], StdI.a)
    StdI.length[1] = print_val_d("Llength", StdI.length[1], StdI.a)
    StdI.direct[0, 0] = print_val_d("Wx", StdI.direct[0, 0], StdI.length[0])
    StdI.direct[0, 1] = print_val_d("Wy", StdI.direct[0, 1], 0.0)
    StdI.direct[1, 0] = print_val_d("Lx", StdI.direct[1, 0], 0.0)
    StdI.direct[1, 1] = print_val_d("Ly", StdI.direct[1, 1], StdI.length[1])

    StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
    StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)

    init_site(StdI, 2)
    StdI.tau[0, 0] = 0.0
    StdI.tau[0, 1] = 0.0
    StdI.tau[0, 2] = 0.0

    # (2) check & store parameters of Hamiltonian
    logger.info("\n  @ Hamiltonian \n")
    not_used_j("J2", StdI.J2All, StdI.J2)
    not_used_j("J2'", StdI.J2pAll, StdI.J2p)
    not_used_d("t2", StdI.t2)
    not_used_d("t2'", StdI.t2p)
    not_used_d("V2", StdI.V2)
    not_used_d("V2'", StdI.V2p)
    not_used_d("K", StdI.K)
    StdI.h = print_val_d("h", StdI.h, 0.0)
    StdI.Gamma = print_val_d("Gamma", StdI.Gamma, 0.0)
    StdI.Gamma_y = print_val_d("Gamma_y", StdI.Gamma_y, 0.0)

    if StdI.model == ModelType.SPIN:
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
        not_used_d("t", StdI.t)
        not_used_d("t0", StdI.t0)
        not_used_d("t1", StdI.t1)
        not_used_d("t'", StdI.tp)
        not_used_d("t0'", StdI.t0p)
        not_used_d("t1'", StdI.t1p)
        not_used_d("t''", StdI.tpp)
        not_used_d("t0''", StdI.t0pp)
        not_used_d("t1''", StdI.t1pp)
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

        if StdI.model == ModelType.HUBBARD:
            not_used_i("2S", StdI.S2)
            not_used_j("J", StdI.JAll, StdI.J)
        else:
            StdI.S2 = print_val_i("2S", StdI.S2, 1)
            input_spin(StdI.J, StdI.JAll, "J")

    logger.info("\n  @ Numerical conditions\n")

    # (3) Set local spin flag and number of sites
    set_local_spin_flags(StdI, StdI.NsiteUC * StdI.NCell)

    # (4) Compute upper limit of Transfer & Interaction and malloc them
    #     nn=2, nnn=2, nnnn=2 → n_bonds=6
    ntransMax = compute_max_interactions(StdI, n_bonds=2 + 2 + 2)
    malloc_interactions(StdI, ntransMax)

    # (5) Set Transfer & Interaction
    # Relative bond table: (dW, dL, site_i, site_j, nn_level, J, t, V)
    _BONDS = (
        (1, 0, 0, 0, 1, StdI.J0, StdI.t0, StdI.V0),       # nn along W
        (0, 1, 0, 0, 1, StdI.J1, StdI.t1, StdI.V1),       # nn along L
        (1, 1, 0, 0, 2, StdI.J0p, StdI.t0p, StdI.V0p),    # nnn W+L
        (1, -1, 0, 0, 2, StdI.J1p, StdI.t1p, StdI.V1p),   # nnn W-L
        (2, 0, 0, 0, 3, StdI.J0pp, StdI.t0pp, StdI.V0pp), # nnnn 2W
        (0, 2, 0, 0, 3, StdI.J1pp, StdI.t1pp, StdI.V1pp), # nnnn 2L
    )
    expand_bonds_2d(StdI, buf, _BONDS)
    return buf.build(StdI) if buf else None


# ---------------------------------------------------------------------------
#  Lattice plugin registration
# ---------------------------------------------------------------------------

from . import LatticePlugin, register_lattice


class SquarePlugin(LatticePlugin):
    """Plugin for the 2D tetragonal (square) lattice."""

    @property
    def name(self) -> str:
        return "tetragonal"

    @property
    def aliases(self) -> list[str]:
        return ["tetragonal", "tetragonallattice", "square", "squarelattice"]

    @property
    def ndim(self) -> int:
        return 2

    def setup(self, StdI: StdIntList) -> "GnuplotData | None":
        """Delegate to the tetragonal() function."""
        return tetragonal(StdI)


register_lattice(SquarePlugin())
