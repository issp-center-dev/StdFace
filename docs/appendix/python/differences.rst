Differences from the C Version
================================

.. note::

   This section describes the **Python reimplementation** of StdFace
   (``python/stdface/``).  The main manual covers the C implementation.

The Python implementation produces output that matches the C version for all
supported solvers and lattices.  Exact text files (e.g. ``namelist.def``) are
byte-identical; files containing floating-point numbers are compared with a
relative tolerance of 1 × 10\ :sup:`-8` (``NaN``/``-NaN`` differences are
also tolerated).  The differences below are limited to invocation and
extensibility.

Comment Character
-----------------

Both versions use ``//`` to introduce comments.  Inline comments are not
supported: everything after ``=`` is taken as the value, so
``t = 1.0 // comment`` would produce a parse error.

.. code-block:: text

   // this line is skipped
   model = Hubbard

Solver Selection
----------------

The C version ships as four separate executables, one per solver:

.. code-block:: bash

   ./hphi_dry.out  stan.in
   ./mvmc_dry.out  stan.in
   ./uhf_dry.out   stan.in
   ./hwave_dry.out stan.in

The Python version is a single command with a ``--solver`` flag:

.. code-block:: bash

   stdface stan.in --solver HPhi
   stdface stan.in --solver mVMC
   stdface stan.in --solver UHF
   stdface stan.in --solver HWAVE

The default solver when ``--solver`` is omitted is HPhi.

Plugin Architecture
-------------------

The C version selects solvers and lattices at **compile time** using
preprocessor macros (``-DHPhi``, ``-DmVMC``, etc.).  Adding a new solver or
lattice requires modifying the C source tree.

The Python version uses a **plugin architecture**: solvers and lattices are
self-registering modules discovered at import time.  A new solver or lattice
can be added without touching existing code — only a new module and an entry in
the discovery list are required.

See :doc:`plugin_tutorial` for step-by-step instructions.

Source-to-Module Mapping
------------------------

The table below shows the correspondence between C source files and their
Python counterparts.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - C Source File
     - Python Module
   * - ``dry.c``
     - ``python/__main__.py``
   * - ``StdFace_main.c``
     - ``stdface/core/stdface_main.py``
   * - ``StdFace_vals.h``
     - ``stdface/core/stdface_vals.py``
   * - ``StdFace_ModelUtil.c/h``
     - ``stdface/core/stdface_model_util.py``
   * - ``version.h``
     - ``stdface/version.py``
   * - ``ChainLattice.c``
     - ``stdface/lattice/chain_lattice.py``
   * - ``SquareLattice.c``
     - ``stdface/lattice/square_lattice.py``
   * - ``Ladder.c``
     - ``stdface/lattice/ladder.py``
   * - ``TriangularLattice.c``
     - ``stdface/lattice/triangular_lattice.py``
   * - ``HoneycombLattice.c``
     - ``stdface/lattice/honeycomb_lattice.py``
   * - ``Kagome.c``
     - ``stdface/lattice/kagome.py``
   * - ``Orthorhombic.c``
     - ``stdface/lattice/orthorhombic.py``
   * - ``FCOrtho.c``
     - ``stdface/lattice/fc_ortho.py``
   * - ``Pyrochlore.c``
     - ``stdface/lattice/pyrochlore.py``
   * - ``Wannier90.c``
     - ``stdface/lattice/wannier90.py``
   * - ``export_wannier90.c``
     - ``stdface/writer/export_wannier90.py``
