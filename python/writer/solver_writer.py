"""Solver-specific output writer classes.

This module provides a class hierarchy for writing Expert-mode definition
files for each supported solver.  It replaces the ``if solver == ...``
dispatch block that was previously in ``stdface_main.stdface_main()``.

Classes
-------
SolverWriter
    Abstract base class defining the writer interface.
HPhiWriter
    Writes definition files for the HPhi solver.
MVMCWriter
    Writes definition files for the mVMC solver.
UHFWriter
    Writes definition files for the UHF solver.
HWaveWriter
    Writes definition files for the H-wave solver.

Functions
---------
get_solver_writer
    Factory function returning the correct writer for a solver name.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

from stdface_vals import StdIntList, SolverType, MethodType, NaN_i
from param_check import print_val_i
from .mvmc_variational import (
    generate_orb,
    proj,
    print_jastrow,
)
from .hphi_writer import (
    print_calc_mod,
    print_excitation,
    print_pump,
)
from .mvmc_writer import (
    print_orb,
    print_orb_para,
    print_gutzwiller,
)
from .common_writer import (
    print_loc_spin,
    print_trans,
    print_namelist,
    print_mod_para,
    print_1_green,
    print_2_green,
    check_output_mode,
    check_mod_para,
)
from .interaction_writer import print_interactions
from .export_wannier90 import export_geometry, export_interaction


class SolverWriter(ABC):
    """Abstract base class for solver-specific output writers.

    Each concrete subclass implements :meth:`write` to produce the set of
    Expert-mode definition files required by a particular solver.

    Parameters
    ----------
    name : str
        Human-readable solver name (e.g. ``"HPhi"``).
    """

    def __init__(self, name: str) -> None:
        self.name = name

    @abstractmethod
    def write(self, StdI: StdIntList) -> None:
        """Write all Expert-mode definition files for this solver.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """


class HPhiWriter(SolverWriter):
    """Writer for the HPhi exact-diagonalisation / Lanczos solver.

    Produces: locspn, trans, interactions, modpara, single/pair (excitation),
    teone/tetwo (time-evolution pump), calcmod, greenone, greentwo, namelist.
    """

    def __init__(self) -> None:
        super().__init__(SolverType.HPhi)

    def write(self, StdI: StdIntList) -> None:
        """Write all HPhi definition files.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """
        print_loc_spin(StdI)
        print_trans(StdI)
        print_interactions(StdI)
        check_mod_para(StdI)
        print_mod_para(StdI)
        print_excitation(StdI)
        if StdI.method == MethodType.TIME_EVOLUTION:
            print_pump(StdI)
        print_calc_mod(StdI)
        check_output_mode(StdI)
        print_1_green(StdI)
        print_2_green(StdI)
        print_namelist(StdI)


class MVMCWriter(SolverWriter):
    """Writer for the mVMC variational Monte Carlo solver.

    Produces: locspn, trans, interactions, modpara, orbitalidx,
    orbitalidxpara/gen, gutzwilleridx, greenone, greentwo, namelist,
    plus Jastrow and projection files.
    """

    def __init__(self) -> None:
        super().__init__(SolverType.mVMC)

    def write(self, StdI: StdIntList) -> None:
        """Write all mVMC definition files.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """
        print_loc_spin(StdI)
        print_trans(StdI)
        print_interactions(StdI)
        check_mod_para(StdI)
        print_mod_para(StdI)

        if StdI.lGC == 0 and (StdI.Sz2 == 0 or StdI.Sz2 == NaN_i):
            StdI.ComplexType = print_val_i("ComplexType", StdI.ComplexType, 0)
        else:
            StdI.ComplexType = print_val_i("ComplexType", StdI.ComplexType, 1)

        generate_orb(StdI)
        proj(StdI)
        print_jastrow(StdI)
        if StdI.lGC == 1 or (StdI.Sz2 != 0 and StdI.Sz2 != NaN_i):
            print_orb_para(StdI)
        print_gutzwiller(StdI)
        print_orb(StdI)
        check_output_mode(StdI)
        print_1_green(StdI)
        print_2_green(StdI)
        print_namelist(StdI)


class UHFWriter(SolverWriter):
    """Writer for the UHF (unrestricted Hartree--Fock) solver.

    Produces: locspn, trans, interactions, modpara, greenone, namelist.
    """

    def __init__(self) -> None:
        super().__init__(SolverType.UHF)

    def write(self, StdI: StdIntList) -> None:
        """Write all UHF definition files.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """
        print_loc_spin(StdI)
        print_trans(StdI)
        print_interactions(StdI)
        check_mod_para(StdI)
        print_mod_para(StdI)
        check_output_mode(StdI)
        print_1_green(StdI)
        print_namelist(StdI)


class HWaveWriter(SolverWriter):
    """Writer for the H-wave solver.

    In ``uhfr`` calc-mode, writes: trans, interactions, modpara, greenone.
    Otherwise, exports Wannier90 geometry and interaction files.
    """

    def __init__(self) -> None:
        super().__init__(SolverType.HWAVE)

    def write(self, StdI: StdIntList) -> None:
        """Write all H-wave definition files.

        Parameters
        ----------
        StdI : StdIntList
            The fully-populated parameter structure.
        """
        if StdI.calcmode == "uhfr":
            print_trans(StdI)
            print_interactions(StdI)
            check_mod_para(StdI)
            check_output_mode(StdI)
            print_1_green(StdI)
        else:
            export_geometry(StdI)
            export_interaction(StdI)


# Registry of solver writers
_SOLVER_WRITERS: dict[SolverType, type[SolverWriter]] = {
    SolverType.HPhi: HPhiWriter,
    SolverType.mVMC: MVMCWriter,
    SolverType.UHF: UHFWriter,
    SolverType.HWAVE: HWaveWriter,
}


def get_solver_writer(solver: str) -> SolverWriter:
    """Return the appropriate solver writer instance.

    Parameters
    ----------
    solver : str
        Solver name.  Must be one of ``"HPhi"``, ``"mVMC"``, ``"UHF"``,
        or ``"HWAVE"``.

    Returns
    -------
    SolverWriter
        An instance of the concrete writer for *solver*.

    Raises
    ------
    ValueError
        If *solver* is not a recognised solver name.
    """
    cls = _SOLVER_WRITERS.get(solver)
    if cls is None:
        raise ValueError(
            f"Unknown solver {solver!r}. "
            f"Expected one of: {', '.join(_SOLVER_WRITERS)}"
        )
    return cls()
