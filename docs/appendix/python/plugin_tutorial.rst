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

Part 1: Adding a New Lattice
-----------------------------

Example: Adding a "dice" (T3) Lattice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dice lattice is a 2D lattice with 3 sites per unit cell.

Step 1: Create the lattice module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create ``python/stdface/lattice/dice.py``:

.. code-block:: python

   """
   Standard mode for the dice (T3) lattice.

   License
   -------
   HPhi-mVMC-StdFace - Common input generator
   Copyright (C) 2015 The University of Tokyo
   (GPL v3)
   """
   from __future__ import annotations

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
   from .site_util import init_site, set_label, set_local_spin_flags, lattice_gp


   def dice(StdI: StdIntList) -> None:
       """Set up the Hamiltonian for the dice (T3) lattice."""
       with lattice_gp(StdI) as fp:

           # --- (1) Geometry ---
           StdI.NsiteUC = 3  # 3 sites per unit cell

           print("  @ Lattice Size & Shape\n")

           StdI.a = print_val_d("a", StdI.a, 1.0)
           StdI.length[0] = print_val_d("Wlength", StdI.length[0], StdI.a)
           StdI.length[1] = print_val_d("Llength", StdI.length[1], StdI.a)
           StdI.direct[0, 0] = print_val_d("Wx", StdI.direct[0, 0], StdI.length[0])
           StdI.direct[0, 1] = print_val_d("Wy", StdI.direct[0, 1], 0.0)
           StdI.direct[1, 0] = print_val_d("Lx", StdI.direct[1, 0],
                                             StdI.length[1] * 0.5)
           StdI.direct[1, 1] = print_val_d("Ly", StdI.direct[1, 1],
                                             StdI.length[1] * math.sqrt(3.0) / 2.0)

           StdI.phase[0] = print_val_d("phase0", StdI.phase[0], 0.0)
           StdI.phase[1] = print_val_d("phase1", StdI.phase[1], 0.0)

           init_site(StdI, fp, 2)

           # tau positions for 3 sites in the unit cell
           StdI.tau[0, :] = [0.0, 0.0, 0.0]
           StdI.tau[1, :] = [1.0 / 3.0, 1.0 / 3.0, 0.0]
           StdI.tau[2, :] = [2.0 / 3.0, 2.0 / 3.0, 0.0]

           # --- (2) Hamiltonian parameters ---
           print("\n  @ Hamiltonian \n")
           # ... (validate and set parameters, same pattern as kagome.py)

           # --- (3) Local spin flags ---
           print("\n  @ Numerical conditions\n")
           set_local_spin_flags(StdI, StdI.W * StdI.L)

           # --- (4) Allocate interactions ---
           ntransMax, nintrMax = compute_max_interactions(StdI, nbond=6)
           malloc_interactions(StdI, ntransMax, nintrMax)

           # --- (5) Build interactions ---
           for iW in range(StdI.W):
               for iL in range(StdI.L):
                   isite = iW + iL * StdI.W
                   add_local_terms(StdI, isite, iL)
                   # add_neighbor_interaction(StdI, fp, iW, iL, ...)


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

       def setup(self, StdI: StdIntList) -> None:
           dice(StdI)


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

Step 1: Create the solver package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   python/stdface/solvers/mysolver/
     __init__.py
     _plugin.py
     writer.py       # (optional) solver-specific output functions

Step 2: Implement the plugin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create ``python/stdface/solvers/mysolver/_plugin.py``:

.. code-block:: python

   """MySolver plugin."""
   from __future__ import annotations

   from ...plugin import SolverPlugin, register
   from ...core.stdface_vals import StdIntList, NaN_i, NaN_d
   from ...core.keyword_parser import (
       store_with_check_dup_i, store_with_check_dup_d,
       store_with_check_dup_s,
   )


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
           return [
               ("MaxIter", NaN_i),
               ("Threshold", NaN_d),
           ]

       @property
       def reset_arrays(self) -> list[tuple[str, object]]:
           return []

       def write_solver_specific(self, StdI: StdIntList) -> None:
           _write_mysolver_config(StdI)

       def write_green(self, StdI: StdIntList) -> None:
           from ...writer.common_writer import print_1_green
           print_1_green(StdI)

       def post_lattice(self, StdI: StdIntList) -> None:
           if StdI.MaxIter == NaN_i:
               StdI.MaxIter = 500


   def _write_mysolver_config(StdI: StdIntList) -> None:
       with open("mysolver_config.def", "w") as fp:
           fp.write(f"MaxIter = {StdI.MaxIter}\n")
           fp.write(f"Threshold = {StdI.Threshold}\n")
           if hasattr(StdI, 'OutputFile'):
               fp.write(f"OutputFile = {StdI.OutputFile}\n")


   register(MySolverPlugin())

Step 3: Make it discoverable
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

Step 4: Add tests
~~~~~~~~~~~~~~~~~

Create ``test/unit/test_mysolver_plugin.py``:

.. code-block:: python

   """Unit tests for the MySolver plugin."""
   from __future__ import annotations

   import pytest
   from stdface.plugin import get_plugin, SolverPlugin


   class TestMySolverPlugin:

       def test_registered(self):
           plugin = get_plugin("MySolver")
           assert isinstance(plugin, SolverPlugin)

       def test_name(self):
           assert get_plugin("MySolver").name == "MySolver"

       def test_keyword_table(self):
           plugin = get_plugin("MySolver")
           assert "maxiter" in plugin.keyword_table
           assert "threshold" in plugin.keyword_table

       def test_reset_scalars(self):
           plugin = get_plugin("MySolver")
           names = [name for name, _ in plugin.reset_scalars]
           assert "MaxIter" in names
           assert "Threshold" in names

Step 5: Use it
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
     - Build lattice geometry and interactions
   * - ``boost(StdI)``
     - method
     - No
     - HPhi Boost mode (default: no-op)

Registration: call ``register_lattice(MyPlugin())`` at module level.

SolverPlugin
^^^^^^^^^^^^

Defined in ``python/stdface/plugin.py``.

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
     - Scalar field reset table
   * - ``reset_arrays``
     - ``list[tuple]`` property
     - Yes
     - Array field reset table
   * - ``write(StdI)``
     - method
     - No
     - Template method (calls steps below)
   * - ``write_locspn(StdI)``
     - method
     - No
     - Write ``locspn.def``
   * - ``write_trans(StdI)``
     - method
     - No
     - Write ``trans.def``
   * - ``write_interactions(StdI)``
     - method
     - No
     - Write interaction files
   * - ``check_and_write_modpara(StdI)``
     - method
     - No
     - Write ``modpara.def``
   * - ``write_solver_specific(StdI)``
     - method
     - No
     - Solver-specific files
   * - ``write_green(StdI)``
     - method
     - No
     - Write Green's function files
   * - ``write_namelist(StdI)``
     - method
     - No
     - Write ``namelist.def``
   * - ``post_lattice(StdI)``
     - method
     - No
     - Hook after lattice construction

Registration: call ``register(MyPlugin())`` at module level.

Template Method Flow (SolverPlugin.write)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   write()
     |-- write_locspn()              # locspn.def
     |-- write_trans()               # trans.def
     |-- write_interactions()        # coulombinter.def, etc.
     |-- check_and_write_modpara()   # modpara.def
     |-- write_solver_specific()     # solver-specific files
     |-- check_output_mode()         # validate output settings
     |-- write_green()               # greenone.def, greentwo.def
     |-- write_namelist()            # namelist.def

Override individual steps to customise output.  Override ``write()`` entirely
for a completely different output sequence.

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
     - ``HWAVE``
