Plugin Tutorial
===============

.. note::

   This section describes the **Python reimplementation** of StdFace
   (``python/stdface/``).  The main manual covers the C implementation.

StdFace's Python implementation uses a plugin architecture where solvers and
lattices are self-registering modules.  This tutorial explains how to add a
new solver or lattice by creating a plugin module and wiring it into the
existing discovery mechanism.

.. contents::
   :local:
   :depth: 2

Overview
--------

How the Plugin System Works
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   User input (stan.in)
     |
     v
   stdface_main() ---> get_lattice("chain") ---> ChainPlugin.setup(StdI)
     |                                                |
     |                                                v
     |                                         (builds geometry & interactions)
     |
     +---> get_plugin("HPhi") ---> HPhiPlugin.write(StdI)
                                        |
                                        v
                                  (writes .def files)

1. The input file specifies ``lattice = ...`` and ``model = ...``.
2. ``stdface_main`` looks up the lattice plugin by name and calls ``setup(StdI)``.
3. After lattice construction, the solver plugin's ``write(StdI)`` generates
   output files.

Both registries use **auto-registration**: when a module is imported it creates
a plugin instance and registers it.  No central dispatch table needs to be
edited.

Data Model: ``StdIntList`` and Per-Solver Configs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``StdIntList`` (``core/stdface_vals.py``) is the parameter structure threaded
through the run.  It holds only **solver-independent** data — lattice geometry,
model couplings, Hamiltonian term lists and common selectors — grouped into
sub-objects (``_lattice`` / ``_model`` / ``_terms`` / ``_w90``) reached through
façade properties.

**Solver-specific** parameters live in a per-solver *config* dataclass
(``solvers/<name>/config.py``) registered with ``register_config``.  Setting
``StdI.solver`` attaches the matching config to ``StdI._solver_cfg``; reads and
writes of solver fields (``StdI.method``, ``StdI.NVMCCalMode``, …) are routed to
it by ``StdIntList.__getattr__`` / ``__setattr__``.  Consequently ``core/``
contains no solver-specific fields, keyword tables or reset tables, and a new
solver is fully contained under ``solvers/<name>/``.

Part 1: Adding a New Lattice
-----------------------------

Example: Adding a "dice" (T3) Lattice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dice lattice is a 2D lattice with 3 sites per unit cell.

Step 1: Create the lattice module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create ``python/stdface/lattice/dice.py``:

.. code-block:: python

   """Standard mode for the dice (T3) lattice."""
   from __future__ import annotations

   import logging
   import math

   from ..core.stdface_vals import StdIntList, ModelType
   from ..core.param_check import print_val_d, print_val_i, not_used_d
   from .input_params import input_spin_nn, input_spin, input_hopp, input_coulomb_v
   from .interaction_builder import (
       compute_max_interactions, malloc_interactions,
       add_neighbor_interaction, add_local_terms,
   )
   from .site_util import (
       init_site, set_local_spin_flags, new_gnuplot_buffer, GnuplotData,
   )

   logger = logging.getLogger(__name__)


   def dice(StdI: StdIntList) -> "GnuplotData | None":
       """Set up the Hamiltonian for the dice (T3) lattice."""
       # gnuplot buffer is None unless this lattice emits a lattice.gp preview
       buf = new_gnuplot_buffer(StdI)

       # --- (1) Geometry ---
       StdI.NsiteUC = 3  # 3 sites per unit cell
       logger.info("  @ Lattice Size & Shape\n")

       StdI.a = print_val_d("a", StdI.a, 1.0)
       StdI.length[0] = print_val_d("Wlength", StdI.length[0], StdI.a)
       StdI.length[1] = print_val_d("Llength", StdI.length[1], StdI.a)
       StdI.direct[0, 0] = print_val_d("Wx", StdI.direct[0, 0], StdI.length[0])
       StdI.direct[0, 1] = print_val_d("Wy", StdI.direct[0, 1], 0.0)
       StdI.direct[1, 0] = print_val_d("Lx", StdI.direct[1, 0], StdI.length[1] * 0.5)
       StdI.direct[1, 1] = print_val_d("Ly", StdI.direct[1, 1],
                                       StdI.length[1] * math.sqrt(3.0) / 2.0)

       StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
       StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)

       init_site(StdI, 2)
       StdI.tau[0, :] = [0.0, 0.0, 0.0]
       StdI.tau[1, :] = [1.0 / 3.0, 1.0 / 3.0, 0.0]
       StdI.tau[2, :] = [2.0 / 3.0, 2.0 / 3.0, 0.0]

       # --- (2) Hamiltonian parameters ---
       logger.info("\n  @ Hamiltonian \n")
       # ... validate and store parameters, same pattern as honeycomb.py ...

       # --- (3) Local spin flags ---
       logger.info("\n  @ Numerical conditions\n")
       set_local_spin_flags(StdI, StdI.NsiteUC * StdI.NCell)

       # --- (4) Allocate interactions (n_bonds = neighbour bonds per cell) ---
       ntransMax = compute_max_interactions(StdI, n_bonds=6)
       malloc_interactions(StdI, ntransMax)

       # --- (5) Build interactions ---
       for kCell in range(StdI.NCell):
           cell_w = StdI.Cell[kCell, 0]
           cell_l = StdI.Cell[kCell, 1]
           isite = StdI.NsiteUC * kCell
           for uc_i in range(StdI.NsiteUC):
               add_local_terms(StdI, isite + uc_i, isite + uc_i)
           # for each bond: add_neighbor_interaction(
           #     StdI, buf, cell_w, cell_l, dW, dL, si, sj, nn, J, t, V)

       return buf.build(StdI) if buf else None


   # ---------------------------------------------------------------------------
   #  Lattice plugin registration
   # ---------------------------------------------------------------------------

   from . import LatticePlugin, register_lattice


   class DicePlugin(LatticePlugin):
       """Plugin for the 2D dice (T3) lattice."""

       @property
       def name(self) -> str:
           return "dice"

       @property
       def aliases(self) -> list[str]:
           return ["dice", "t3", "dicelattice"]

       @property
       def ndim(self) -> int:
           return 2

       def setup(self, StdI: StdIntList) -> "GnuplotData | None":
           return dice(StdI)


   register_lattice(DicePlugin())

Key points:

- The ``LatticePlugin`` subclass defines ``name``, ``aliases``, ``ndim``, and
  ``setup()``.
- ``register_lattice()`` is called at module level — the plugin registers itself
  when the module is imported.
- If the lattice supports HPhi Boost mode, override ``boost(self, StdI)``.

Step 2: Register the module for auto-discovery
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit ``python/stdface/lattice/__init__.py`` and add your module to
``_discover_lattices()``:

.. code-block:: python

   def _discover_lattices() -> None:
       """Import all lattice modules to trigger auto-registration."""
       from . import chain_lattice, square_lattice, ladder, triangular_lattice
       from . import honeycomb_lattice, kagome, orthorhombic, fc_ortho, pyrochlore
       from . import wannier90
       from . import dice  # <-- Add this line

Now ``lattice = dice`` (or ``lattice = t3``) will work in input files.

Step 3: Add tests
~~~~~~~~~~~~~~~~~

Create ``test/unit/test_dice.py``:

.. code-block:: python

   """Unit tests for the dice lattice plugin."""
   from __future__ import annotations

   import pytest
   from stdface.lattice import get_lattice, LatticePlugin


   class TestDicePlugin:

       def test_registered(self):
           plugin = get_lattice("dice")
           assert isinstance(plugin, LatticePlugin)

       def test_name(self):
           assert get_lattice("dice").name == "dice"

       def test_ndim(self):
           assert get_lattice("dice").ndim == 2

       def test_aliases(self):
           p = get_lattice("dice")
           assert get_lattice("t3") is p
           assert get_lattice("dicelattice") is p

       def test_setup_callable(self):
           assert callable(get_lattice("dice").setup)

Step 4: Verify
~~~~~~~~~~~~~~

.. code-block:: bash

   # Unit tests
   python3 -m pytest test/unit/ -q

   # Integration test
   /path/to/StdFace/stdface stan.in

Part 2: Adding a New Solver
-----------------------------

Example: Adding a "MySolver" Solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A solver is **self-contained**: every solver-specific parameter lives in a
per-solver *config* dataclass (``config.py``), and the plugin
(``_plugin.py``) owns the input keyword/reset tables and the output.  Nothing
in ``core/`` has to be edited to add a solver.

Step 1: Create the solver package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   python/stdface/solvers/mysolver/
     __init__.py
     config.py       # per-solver field container (a @dataclass)
     _plugin.py      # plugin: keyword/reset tables + output

Step 2: Declare the solver's fields (config)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Solver-specific fields no longer live on ``StdIntList``; each solver owns a
``@dataclass`` that is attached to ``StdI`` for the duration of a run.  Create
``python/stdface/solvers/mysolver/config.py``:

.. code-block:: python

   """MySolver field container."""
   from __future__ import annotations

   from dataclasses import dataclass


   @dataclass
   class MySolverConfig:
       MaxIter: int | None = None
       Threshold: float | None = None
       OutputFile: str | None = None

Unset optional fields use ``None`` as the sentinel (numpy integer arrays use
the ``NaN_i`` sentinel instead, since ``None`` cannot be stored in them).

Step 3: Implement the plugin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create ``python/stdface/solvers/mysolver/_plugin.py``.  Reading or writing
``StdI.MaxIter`` transparently reaches the attached ``MySolverConfig`` (via
``StdIntList.__getattr__`` / ``__setattr__``):

.. code-block:: python

   """MySolver plugin."""
   from __future__ import annotations

   from ...plugin import SolverPlugin, register, register_config
   from ...core.stdface_vals import StdIntList
   from ...core.keyword_parser import (
       store_with_check_dup_i, store_with_check_dup_d,
       store_with_check_dup_s,
   )
   from .config import MySolverConfig


   class MySolverPlugin(SolverPlugin):
       """Plugin for the MySolver backend."""

       @property
       def name(self) -> str:
           return "MySolver"

       @property
       def keyword_table(self) -> dict[str, tuple]:
           return {
               "maxiter":    (store_with_check_dup_i, "MaxIter"),
               "threshold":  (store_with_check_dup_d, "Threshold"),
               "outputfile": (store_with_check_dup_s, "OutputFile"),
           }

       @property
       def reset_scalars(self) -> list[tuple[str, object]]:
           return [("MaxIter", None), ("Threshold", None)]

       @property
       def reset_arrays(self) -> list[tuple[str, object]]:
           return []

       def set_defaults(self, StdI: StdIntList) -> None:
           if StdI.MaxIter is None:
               StdI.MaxIter = 500

       def write(self, StdI: StdIntList) -> None:
           with open("mysolver_config.def", "w") as fp:
               fp.write(f"MaxIter = {StdI.MaxIter}\n")
               fp.write(f"Threshold = {StdI.Threshold}\n")
               if StdI.OutputFile is not None:
                   fp.write(f"OutputFile = {StdI.OutputFile}\n")


   # Auto-register on import: plugin under its name, config under the same name.
   register(MySolverPlugin())
   register_config("MySolver", MySolverConfig)

.. note::

   Solvers that emit the standard HPhi/mVMC-style files (``locspn.def`` /
   ``trans.def`` / ``modpara.def`` / ``namelist.def`` / Green's functions)
   should subclass :class:`ExpertModeSolverPlugin` instead of
   ``SolverPlugin`` and implement :meth:`modpara_lines` (and optionally
   ``namelist_entries`` / ``has_two_body_green`` / ``write_solver_files``);
   the base class assembles them into an ``ExpertModeOutput`` and writes it.
   See ``solvers/uhf/`` for a minimal example.

Step 4: Make it discoverable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create ``python/stdface/solvers/mysolver/__init__.py``:

.. code-block:: python

   """MySolver plugin package."""
   from . import _plugin  # noqa: F401 — triggers auto-registration

Edit ``python/stdface/solvers/__init__.py``:

.. code-block:: python

   """Built-in solver plugins for StdFace."""
   from __future__ import annotations

   from . import hphi, mvmc, uhf, hwave  # noqa: F401
   from . import mysolver  # noqa: F401  <-- Add this line

Step 5: Add tests
~~~~~~~~~~~~~~~~~

Create ``test/unit/test_mysolver_plugin.py``:

.. code-block:: python

   """Unit tests for the MySolver plugin."""
   from __future__ import annotations

   import pytest
   from stdface.plugin import get_plugin, SolverPlugin
   from stdface.core.stdface_vals import StdIntList


   class TestMySolverPlugin:

       def test_registered(self):
           assert isinstance(get_plugin("MySolver"), SolverPlugin)

       def test_keyword_table(self):
           kw = get_plugin("MySolver").keyword_table
           assert "maxiter" in kw and "threshold" in kw

       def test_config_attached_by_solver(self):
           # Setting the solver auto-attaches MySolverConfig; its fields resolve.
           s = StdIntList()
           s.solver = "MySolver"
           s.MaxIter = 10
           assert s.MaxIter == 10

Step 6: Use it
~~~~~~~~~~~~~~

.. code-block:: bash

   ./stdface stan.in --solver MySolver

Plugin API Reference
--------------------

LatticePlugin
^^^^^^^^^^^^^

Defined in ``python/stdface/lattice/__init__.py``.

.. list-table::
   :header-rows: 1
   :widths: 20 15 10 55

   * - Member
     - Type
     - Required
     - Description
   * - ``name``
     - ``str`` property
     - Yes
     - Canonical name (e.g. ``"chain"``)
   * - ``aliases``
     - ``list[str]`` property
     - Yes
     - All recognized aliases
   * - ``ndim``
     - ``int`` property
     - Yes
     - Spatial dimensions (1, 2, or 3)
   * - ``setup(StdI)``
     - method
     - Yes
     - Build geometry/interactions; return ``GnuplotData | None``
   * - ``boost(StdI)``
     - method
     - No
     - HPhi Boost mode (default: no-op)

Registration: call ``register_lattice(MyPlugin())`` at module level.

SolverPlugin
^^^^^^^^^^^^

Defined in ``python/stdface/plugin.py``.  A solver implements the four
required members plus ``write``; the optional hooks default to no-ops.

.. list-table::
   :header-rows: 1
   :widths: 30 15 10 45

   * - Member
     - Type
     - Required
     - Description
   * - ``name``
     - ``str`` property
     - Yes
     - Canonical name (e.g. ``"HPhi"``)
   * - ``keyword_table``
     - ``dict`` property
     - Yes
     - Input keyword dispatch table
   * - ``reset_scalars``
     - ``list[tuple]`` property
     - Yes
     - Scalar field reset table (``(name, value)``)
   * - ``reset_arrays``
     - ``list[tuple]`` property
     - Yes
     - Array-fill reset table (``(name, fill_value)``)
   * - ``write(StdI)``
     - method
     - Yes
     - Emit this solver's output files
   * - ``set_defaults(StdI)``
     - method
     - No
     - Apply default parameter values (called from ``check_mod_para``)
   * - ``post_lattice(StdI)``
     - method
     - No
     - Hook after lattice construction (e.g. HPhi Boost / LargeValue)
   * - ``init_fields(StdI)``
     - method
     - No
     - Hook to attach extra dynamic attributes
   * - ``validate(StdI)``
     - method
     - No
     - Reject unsupported parameter combinations (raise ``ValueError``)

Registration: call ``register(MyPlugin())`` and ``register_config("Name",
MyConfig)`` at module level.

ExpertModeSolverPlugin
^^^^^^^^^^^^^^^^^^^^^^^

Base class for solvers that emit the HPhi/mVMC-style files (HPhi, mVMC, UHF).
It implements ``write`` via ``build_output`` and adds:

.. list-table::
   :header-rows: 1
   :widths: 30 15 10 45

   * - Member
     - Type
     - Required
     - Description
   * - ``modpara_lines(StdI)``
     - method
     - Yes
     - Solver-specific body lines of ``modpara.def``
   * - ``namelist_entries(StdI)``
     - method
     - No
     - Extra ``(keyword, filename)`` namelist entries (default: none)
   * - ``has_two_body_green(StdI)``
     - method
     - No
     - List ``greentwo.def`` in the namelist (default: ``True``)
   * - ``write_solver_files(StdI)``
     - method
     - No
     - Emit solver-specific extra files (e.g. mVMC variational files)
   * - ``build_output(StdI)``
     - method
     - No
     - Assemble an ``ExpertModeOutput`` (override for a custom container)

Output Flow (ExpertModeSolverPlugin)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Output is built into a data container first, then written — there is no
per-file template method.  ``write()`` calls ``build_output()``:

.. code-block:: text

   build_output()                     # -> ExpertModeOutput
     |-- check_mod_para() / set_defaults()
     |-- build_loc_spn()              # locspn.def
     |-- build_trans()                # trans.def
     |-- build_interactions()         # coulomb*/hund/exchange/... (sets L* flags)
     |-- build_modpara()              # modpara.def
     |-- build_green_one() / build_green_two()
     |-- write_solver_files()         # solver-specific extra files
     |-- build_namelist()             # namelist.def
   write() = build_output(StdI).write()

A plain ``SolverPlugin`` (no modpara/namelist) just overrides ``write``
directly, or builds its own output container.  See ``solvers/hwave/`` for a
plugin that returns either an ``ExpertModeOutput`` (UHFR ``.def``) or a
``WannierModeOutput`` (UHFK Wannier90) depending on ``calcmode``.

Existing Plugins Reference
--------------------------

Lattice Plugins
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 30 15 20 10

   * - Plugin Class
     - Module
     - name
     - aliases
     - ndim
   * - ``ChainPlugin``
     - ``chain_lattice.py``
     - ``chain``
     - chain, chainlattice
     - 1
   * - ``LadderPlugin``
     - ``ladder.py``
     - ``ladder``
     - ladder, ladderlattice
     - 1
   * - ``SquarePlugin``
     - ``square_lattice.py``
     - ``tetragonal``
     - tetragonal, square, squarelattice
     - 2
   * - ``TriangularPlugin``
     - ``triangular_lattice.py``
     - ``triangular``
     - triangular, triangularlattice
     - 2
   * - ``HoneycombPlugin``
     - ``honeycomb_lattice.py``
     - ``honeycomb``
     - honeycomb, honeycomblattice
     - 2
   * - ``KagomePlugin``
     - ``kagome.py``
     - ``kagome``
     - kagome, kagomelattice
     - 2
   * - ``OrthorhombicPlugin``
     - ``orthorhombic.py``
     - ``orthorhombic``
     - orthorhombic, cubic, simplecubic
     - 3
   * - ``FCOrthoPlugin``
     - ``fc_ortho.py``
     - ``fco``
     - fco, fccubic, fcc
     - 3
   * - ``PyrochlorePlugin``
     - ``pyrochlore.py``
     - ``pyrochlore``
     - pyrochlore
     - 3
   * - ``Wannier90Plugin``
     - ``wannier90.py``
     - ``wannier90``
     - wannier90
     - 3

Lattices with Boost support: chain, honeycomb, kagome, ladder.

Solver Plugins
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Plugin Class
     - Module
     - name
   * - ``HPhiPlugin``
     - ``solvers/hphi/_plugin.py``
     - ``HPhi``
   * - ``MVMCPlugin``
     - ``solvers/mvmc/_plugin.py``
     - ``mVMC``
   * - ``UHFPlugin``
     - ``solvers/uhf/_plugin.py``
     - ``UHF``
   * - ``HWavePlugin``
     - ``solvers/hwave/_plugin.py``
     - ``HWAVE`` (also ``UHFR`` / ``UHFK``)

Each solver also has a ``config.py`` (``HPhiConfig`` / ``MVMCConfig`` /
``UHFConfig`` / ``HWaveConfig``).  ``HWavePlugin`` is registered under
``HWAVE`` plus the ``UHFR`` / ``UHFK`` aliases and picks its output mode
(real-space ``.def`` vs Wannier90) from ``calcmode`` at write time.
