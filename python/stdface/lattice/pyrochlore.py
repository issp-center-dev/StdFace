"""
Standard mode for the pyrochlore lattice.

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
    mag_field, general_j, hubbard_local,
    expand_bonds_3d,
)
from .site_util import init_site, set_local_spin_flags


logger = logging.getLogger(__name__)


def pyrochlore(StdI: StdIntList) -> None:
    """Setup a Hamiltonian for the pyrochlore lattice.

    Parameters
    ----------
    StdI : StdIntList
        Structure containing model parameters and lattice information.
        Modified in-place.
    """
    # (1) Compute the shape of the super-cell and sites in the super-cell
    StdI.NsiteUC = 4

    logger.info("  @ Lattice Size & Shape\n")

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

    init_site(StdI, 3)
    StdI.tau[0, 0] = 0.0
    StdI.tau[0, 1] = 0.0
    StdI.tau[0, 2] = 0.0
    StdI.tau[1, 0] = 0.5
    StdI.tau[1, 1] = 0.0
    StdI.tau[1, 2] = 0.0
    StdI.tau[2, 0] = 0.0
    StdI.tau[2, 1] = 0.5
    StdI.tau[2, 2] = 0.0
    StdI.tau[3, 0] = 0.0
    StdI.tau[3, 1] = 0.0
    StdI.tau[3, 2] = 0.5

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
        input_spin_nn(StdI.J, StdI.JAll, StdI.J0p, StdI.J0pAll, "J0'")
        input_spin_nn(StdI.J, StdI.JAll, StdI.J1p, StdI.J1pAll, "J1'")
        input_spin_nn(StdI.J, StdI.JAll, StdI.J2p, StdI.J2pAll, "J2'")

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
        not_used_d("t''", StdI.tpp)
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
        StdI.t0p = input_hopp(StdI.t, StdI.t0p, "t0'")
        StdI.t1p = input_hopp(StdI.t, StdI.t1p, "t1'")
        StdI.t2p = input_hopp(StdI.t, StdI.t2p, "t2'")
        StdI.V0 = input_coulomb_v(StdI.V, StdI.V0, "V0")
        StdI.V1 = input_coulomb_v(StdI.V, StdI.V1, "V1")
        StdI.V2 = input_coulomb_v(StdI.V, StdI.V2, "V2")
        StdI.V0p = input_coulomb_v(StdI.V, StdI.V0p, "V0'")
        StdI.V1p = input_coulomb_v(StdI.V, StdI.V1p, "V1'")
        StdI.V2p = input_coulomb_v(StdI.V, StdI.V2p, "V2'")

        not_used_j("J0", StdI.J0All, StdI.J0)
        not_used_j("J1", StdI.J1All, StdI.J1)
        not_used_j("J2", StdI.J2All, StdI.J2)
        not_used_j("J0'", StdI.J0pAll, StdI.J0p)
        not_used_j("J1'", StdI.J1pAll, StdI.J1p)
        not_used_j("J2'", StdI.J2pAll, StdI.J2p)
        not_used_j("J''", StdI.JppAll, StdI.Jpp)
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
    #     nn=12 → n_bonds=12
    ntransMax = compute_max_interactions(StdI, n_bonds=12)
    malloc_interactions(StdI, ntransMax)

    # (5) Set Transfer & Interaction
    def _local(StdI, kCell):
        isite = StdI.NsiteUC * kCell
        if StdI.model == ModelType.KONDO:
            isite += StdI.nsite // 2

        if StdI.model == ModelType.SPIN:
            for uc_i in range(StdI.NsiteUC):
                mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, isite + uc_i)
                general_j(StdI, StdI.D, StdI.S2, StdI.S2, isite + uc_i, isite + uc_i)
        else:
            for uc_i in range(StdI.NsiteUC):
                hubbard_local(StdI, StdI.mu, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, StdI.U, isite + uc_i)
            if StdI.model == ModelType.KONDO:
                jsite = StdI.NsiteUC * kCell
                for uc_i in range(StdI.NsiteUC):
                    general_j(StdI, StdI.J, 1, StdI.S2, isite + 3, jsite + uc_i)
                    mag_field(StdI, StdI.S2, -StdI.h, -StdI.Gamma, -StdI.Gamma_y, jsite + uc_i)

    # Relative bond table: (dW, dL, dH, site_i, site_j, J, t, V)
    _BONDS = (
        # Intra-cell
        (0, 0, 0, 0, 1, StdI.J0, StdI.t0, StdI.V0),    # along W
        (0, 0, 0, 0, 2, StdI.J1, StdI.t1, StdI.V1),    # along L
        (0, 0, 0, 0, 3, StdI.J2, StdI.t2, StdI.V2),    # along H
        (0, 0, 0, 2, 3, StdI.J0p, StdI.t0p, StdI.V0p), # along L-H
        (0, 0, 0, 3, 1, StdI.J1p, StdI.t1p, StdI.V1p), # along H-W
        (0, 0, 0, 1, 2, StdI.J2p, StdI.t2p, StdI.V2p), # along W-L
        # Inter-cell
        (1, 0, 0, 1, 0, StdI.J0, StdI.t0, StdI.V0),      # along W
        (0, 1, 0, 2, 0, StdI.J1, StdI.t1, StdI.V1),      # along L
        (0, 0, 1, 3, 0, StdI.J2, StdI.t2, StdI.V2),      # along H
        (0, -1, 1, 3, 2, StdI.J0p, StdI.t0p, StdI.V0p),  # along L-H
        (1, 0, -1, 1, 3, StdI.J1p, StdI.t1p, StdI.V1p),  # along H-W
        (-1, 1, 0, 2, 1, StdI.J2p, StdI.t2p, StdI.V2p),  # along W-L
    )
    expand_bonds_3d(StdI, _BONDS, _local)


# ---------------------------------------------------------------------------
#  Lattice plugin registration
# ---------------------------------------------------------------------------

from . import LatticePlugin, register_lattice

_pyrochlore_setup = pyrochlore


class PyrochlorePlugin(LatticePlugin):
    """Plugin for the 3D pyrochlore lattice."""

    @property
    def name(self) -> str:
        return "pyrochlore"

    @property
    def aliases(self) -> list[str]:
        return ["pyrochlore"]

    @property
    def ndim(self) -> int:
        return 3

    def setup(self, StdI: StdIntList) -> None:
        """Delegate to the pyrochlore() function."""
        _pyrochlore_setup(StdI)


register_lattice(PyrochlorePlugin())
