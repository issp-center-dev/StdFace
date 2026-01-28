"""Output file generation subpackage.

This package contains modules for writing Expert-mode definition files
(``.def``) for all supported solvers (HPhi, mVMC, UHF, H-wave).
"""
from .solver_writer import get_solver_writer
from .interaction_writer import print_interactions
