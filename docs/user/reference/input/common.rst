Common Parameters
=================

.. contents::
   :local:
   :depth: 2

Input File Format
-----------------

Basic Syntax
^^^^^^^^^^^^

StdFace input files use a simple ``key = value`` format:

.. code-block:: text

   model = "Hubbard"
   lattice = square
   W = 2
   L = 2
   t = 1.0
   U = 1.0

Formatting rules:

- One key-value pair per line
- Keys and values separated by ``=``
- Whitespace around ``=`` is optional (both ``W=2`` and ``W = 2`` are accepted)
- String values may be quoted (``"Hubbard"``) or unquoted (``square``)
- Numeric values are written without quotes

Comment Syntax
^^^^^^^^^^^^^^

Lines beginning with ``//`` are treated as comments and ignored:

.. code-block:: text

   // This is a comment line
   model = "Hubbard"
   lattice = square

Empty lines are also skipped.  **Inline comments are not supported** — the
parser splits each line on the first ``=`` sign, so a trailing ``// comment``
would be included in the value and cause a parse error.

Ordering
^^^^^^^^

Key ordering within the input file is not required.  Parameters may appear
in any order.

Parameters
----------

The following parameters are accepted by all solvers.

.. list-table::
   :header-rows: 1
   :widths: 15 20 65

   * - Key
     - Example
     - Description
   * - ``model``
     - ``"Hubbard"``
     - Physical model type
   * - ``lattice``
     - ``square``
     - Lattice geometry (see :doc:`../lattice/index`)
   * - ``W``
     - ``2``
     - Lattice dimension in W direction
   * - ``L``
     - ``2``
     - Lattice dimension in L direction
   * - ``t``
     - ``1.0``
     - Nearest-neighbor hopping amplitude
   * - ``U``
     - ``1.0``
     - On-site Coulomb interaction strength

``model``
^^^^^^^^^

Supported values: ``"Hubbard"``, ``"Spin"``, ``"Kondo"``.
Grand-canonical variants (``"HubbardGC"``, ``"SpinGC"``, ``"KondoGC"``) set
the grand-canonical ensemble flag (``lGC = 1``).

``lattice``
^^^^^^^^^^^

Specifies the lattice geometry.  Supported values (canonical names and
aliases):

- ``chain`` / ``chainlattice``
- ``ladder`` / ``ladderlattice``
- ``square`` / ``squarelattice`` / ``tetragonal``
- ``triangular`` / ``triangularlattice``
- ``honeycomb`` / ``honeycomblattice``
- ``kagome`` / ``kagomelattice``
- ``orthorhombic`` / ``simplecubic``
- ``fco`` / ``fcc``
- ``pyrochlore``
- ``wannier90``

See :doc:`../lattice/index` for per-lattice details.

``t``
^^^^^

Nearest-neighbor hopping amplitude.  Not used when
``lattice = "wannier90"`` (hopping is read from ``zvo_hr.dat``).

``U``
^^^^^

On-site Coulomb interaction strength.  Not used when
``lattice = "wannier90"`` (interaction is read from ``zvo_ur.dat``).

Cross-Solver Parameters
-----------------------

The following parameters appear in all solver modes but behave differently
depending on the solver.

``nelec`` — number of electrons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Solver
     - Behaviour
   * - HPhi
     - Required for Hubbard and Kondo models in canonical ensemble
       (``lGC = 0``, the default).  Not needed for spin models.
   * - mVMC
     - Always required.
   * - UHF
     - Always required.
   * - H-wave
     - Always required.

``2Sz`` — total spin projection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Specifies 2×Sz (twice the z-component of total spin).  For example,
``2Sz = 0`` is the Sz = 0 sector; ``2Sz = 1`` is Sz = 1/2.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Solver
     - Behaviour
   * - HPhi
     - Required for spin models.  Optional for Hubbard/Kondo models (defaults
       to 0).  Ignored when grand-canonical ensemble (``lGC = 1``) is used.
   * - mVMC
     - Optional.  When set (and non-zero), ``orbitalidxpara.def`` is generated
       instead of ``orbitalidx.def``.
   * - UHF
     - Optional.
   * - H-wave
     - Optional.
