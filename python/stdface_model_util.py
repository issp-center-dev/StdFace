"""
Utility functions for constructing lattice models in Standard mode.

This module was formerly a monolithic utility library (~1,700 lines).
It has been decomposed into focused modules during Phase 2 refactoring.
All public symbols are re-exported here for backward compatibility so
that existing lattice modules can continue to use::

    from stdface_model_util import init_site, find_site, ...

without modification.

Re-exported modules
-------------------
param_check
    Parameter validation and printing utilities.
input_params
    Input parameter resolution helpers.
geometry_output
    Geometry and structure output functions.
interaction_builder
    Hamiltonian term builder functions and array allocation.
site_util
    Super-cell initialisation, site folding, finding, and labelling.

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

from param_check import (  # noqa: F401 – re-exported for backward compatibility
    exit_program,
    print_val_d,
    print_val_dd,
    print_val_c,
    print_val_i,
    not_used_d,
    not_used_j,
    not_used_i,
    required_val_i,
)
from lattice.input_params import (  # noqa: F401 – re-exported for backward compatibility
    input_spin_nn,
    input_spin,
    input_coulomb_v,
    input_hopp,
)
from lattice.geometry_output import (  # noqa: F401 – re-exported for backward compatibility
    print_xsf,
    print_geometry,
)
from lattice.interaction_builder import (  # noqa: F401 – re-exported for backward compatibility
    trans,
    hopping,
    hubbard_local,
    mag_field,
    intr,
    general_j,
    coulomb,
    malloc_interactions,
)
from lattice.site_util import (  # noqa: F401 – re-exported for backward compatibility
    _fold_site,
    init_site,
    find_site,
    set_label,
)
