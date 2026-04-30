Input / Output Reference
========================

This document describes the input file format for StdFace and the output files
it generates.

Overview
--------

StdFace is an input file generator for quantum lattice model solvers. It reads
a configuration file specifying the physical model and lattice geometry, then
generates solver-specific input files.

**Scope**:

- **Input**: StdFace reads a single text configuration file (commonly named
  ``stan.in``) containing model and lattice parameters.
- **Output**: StdFace writes one or more input files for the target solver
  (HPhi, mVMC, UHF, or H-wave) in the current working directory.

The specific output files depend on the solver mode, which is determined at
compile time.  See `Output Files`_ for a complete list.

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

Observed formatting rules (from sample files):

- One key-value pair per line
- Keys and values separated by ``=``
- Whitespace around ``=`` is optional (both ``W=2`` and ``W = 2`` appear in samples)
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

Key ordering within the input file is unspecified. The sample files do not
suggest any required ordering.

Common Parameters
-----------------

The following keys are shared across solvers.

===========  =============  ================================================
Key          Example        Description
===========  =============  ================================================
``model``    ``"Hubbard"``  Physical model type
``lattice``  ``square``     Lattice geometry
``W``        ``2``          Lattice dimension in W direction
``L``        ``2``          Lattice dimension in L direction
``t``        ``1.0``        Nearest-neighbor hopping amplitude
``U``        ``1.0``        On-site Coulomb interaction strength
===========  =============  ================================================

``model``
   Physical model type.  Supported values include ``"Hubbard"``, ``"Spin"``,
   and ``"Kondo"``.

``lattice``
   Lattice geometry.  Supported values:

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

.. _lattice-types:

Lattice Types
^^^^^^^^^^^^^

==================  ===========  ============================================
Lattice             Dimensions   Description
==================  ===========  ============================================
``chain``           1D           Linear chain
``ladder``          1D           Two-leg or multi-leg ladder
``square``          2D           Square lattice
``triangular``      2D           Triangular lattice
``honeycomb``       2D           Honeycomb (hexagonal) lattice
``kagome``          2D           Kagome lattice
``tetragonal``      3D           Tetragonal lattice
``orthorhombic``    3D           Orthorhombic/cubic lattice
``fcortho``         3D           Face-centered orthorhombic lattice
``pyrochlore``      3D           Pyrochlore lattice
``wannier90``       Any          Import from Wannier90 format files
==================  ===========  ============================================

See the developer documentation for implementation details of each lattice
constructor.

``W``, ``L``
   Lattice dimensions.  Interpretation depends on lattice type.

``t``
   Nearest-neighbor hopping amplitude.  Not used when
   ``lattice = "wannier90"`` (hopping comes from ``zvo_hr.dat``).

``U``
   On-site Coulomb interaction strength.  Not used when
   ``lattice = "wannier90"`` (interaction comes from ``zvo_ur.dat``).

Solver-Specific Parameters
--------------------------

The following parameters are interpreted differently — or only apply — for
particular solver modes.

.. _method-hphi:

``method`` — HPhi only
^^^^^^^^^^^^^^^^^^^^^^^

Specifies the diagonalisation algorithm.  Ignored by mVMC, UHF, and H-wave.

====================  =========================================================
Value                 Algorithm
====================  =========================================================
``"lanczos"``         Lanczos method
``"lanczosenergy"``   Lanczos with eigenvector output
``"tpq"``             Thermal Pure Quantum states
``"fulldiag"``        Full (exact) diagonalisation (aliases: ``"direct"``,
                      ``"alldiag"``)
``"cg"``              Conjugate Gradient
``"timeevolution"``   Real-time evolution (aliases: ``"te"``,
                      ``"time-evolution"``)
``"ctpq"``            Canonical TPQ
====================  =========================================================

``exct`` — HPhi only
^^^^^^^^^^^^^^^^^^^^^

Number of eigenvalues/eigenvectors to compute (ground state + excited states).
Default: ``1`` (ground state only).  Ignored by mVMC, UHF, and H-wave.

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

HPhi Boost Mode
^^^^^^^^^^^^^^^

Boost mode is a high-performance algorithm for S=1/2 spin models that
exploits a 6-spin interaction structure to accelerate the Lanczos
diagonalisation.  Internally, HPhi reorganises the Hilbert space around
*pivot sites* and applies quantum-number projection, achieving better cache
efficiency than the standard algorithm on supported lattices.

**Enabling Boost mode**

Set ``model`` to one of:

- ``"spingcboost"`` — grand-canonical spin with Boost
- ``"spingccma"`` — grand-canonical spin with CMA variant

.. code-block:: text

   model   = "spingcboost"
   lattice = chain
   L       = 16

**Supported lattices and constraints**

+---------------+--------------------------------------------------+
| Lattice       | Constraints                                      |
+===============+==================================================+
| ``chain``     | ``S2 = 1`` (S=1/2); ``L`` must be a multiple of 8|
+---------------+--------------------------------------------------+
| ``ladder``    | ``S2 = 1``; ``W = 2``; ``L`` even and ≥ 4        |
+---------------+--------------------------------------------------+
| ``honeycomb`` | ``S2 = 1``                                       |
+---------------+--------------------------------------------------+
| ``kagome``    | ``S2 = 1``                                       |
+---------------+--------------------------------------------------+

Other lattices do not support Boost mode and will produce an error.

**Additional output file**

When Boost mode is active, ``boost.def`` is generated in addition to the
standard HPhi files.  It contains the magnetic field values, exchange
interaction matrices, and the pivot-site topology (pivot count, 6-spin
pair lists, spin-shift information) required by HPhi's Boost kernel.

Mode-Specific Inputs
--------------------

Wannier90-Based Input
^^^^^^^^^^^^^^^^^^^^^

When ``lattice = "wannier90"`` is specified, StdFace reads lattice and
interaction data from external files in Wannier90 format.

**Observed in**: ``samples/hubbard/wannier/``

Required Files
""""""""""""""

The wannier sample directory contains:

1. ``stan.in`` - Main configuration file
2. ``zvo_geom.dat`` - Geometry data (lattice vectors and atomic positions)
3. ``zvo_hr.dat`` - Hopping (transfer) integrals
4. ``zvo_ur.dat`` - On-site interaction data

Main Input (stan.in)
""""""""""""""""""""

.. code-block:: text

   model = "Hubbard"
   lattice = "wannier90"
   W=2
   L=2
   method = "CG"
   2Sz = 0
   nelec = 4
   exct = 1

Note: ``t`` and ``U`` keys are absent; these values come from the external
Wannier90 files.

Geometry File (zvo_geom.dat)
""""""""""""""""""""""""""""

Format (observed):

.. code-block:: text

     1.0000000000   0.0000000000   0.0000000000
     0.0000000000   1.0000000000   0.0000000000
     0.0000000000   0.0000000000   1.0000000000
     1
     0.5000000000 0.5000000000 0.5000000000

- Lines 1-3: Lattice vectors (3x3 matrix)
- Line 4: Number of atoms in unit cell
- Lines 5+: Fractional coordinates of atoms

Hopping File (zvo_hr.dat)
"""""""""""""""""""""""""

Format (observed):

.. code-block:: text

   wannier90 format for vmcdry.out or HPhi -sdry
            1
            9
       1    1    1    1    1    1    1    1    1
      -1   -1    0    1    1    0.0  0.0
      -1    0    0    1    1   -1.0  0.0
      ...

- Line 1: Header comment
- Line 2: Number of Wannier functions
- Line 3: Number of R vectors
- Line 4: Degeneracy weights
- Lines 5+: (R_x, R_y, R_z, orbital_i, orbital_j, Re(t), Im(t))

Interaction File (zvo_ur.dat)
"""""""""""""""""""""""""""""

Same format as zvo_hr.dat. The entry at R=(0,0,0) with value 1.0 represents
the on-site Hubbard U interaction.

File Relationships
""""""""""""""""""

.. code-block:: text

   stan.in (lattice = "wannier90")
       |
       +-- zvo_geom.dat  (lattice geometry)
       +-- zvo_hr.dat    (hopping parameters)
       +-- zvo_ur.dat    (interaction parameters)

All files must be present in the same directory when running StdFace with
``lattice = "wannier90"``.

Output Files
------------

StdFace writes input files for the target solver into the current working
directory.  The solver mode is fixed at compile time:

==============  ===================  =====================================
Build Option    Executable           Target Solver
==============  ===================  =====================================
``-DHPHI=ON``   ``hphi_dry.out``     HPhi (Exact Diagonalisation)
``-DMVMC=ON``   ``mvmc_dry.out``     mVMC (Variational Monte Carlo)
``-DUHF=ON``    ``uhf_dry.out``      UHF (Unrestricted Hartree-Fock)
``-DHWAVE=ON``  ``hwave_dry.out``    H-wave
==============  ===================  =====================================

Files Common to All Solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following files are always generated regardless of solver:

==================  =========================================================
File                Contents
==================  =========================================================
``namelist.def``    List of all generated definition files
``modpara.def``     Model parameters (system size, electron count, etc.)
``locspn.def``      Local spin configuration for each site
``trans.def``       One-body transfer terms
``lattice.gp``      Gnuplot script for visualising the lattice
``lattice.xsf``     XCrySDen structure file for the lattice
==================  =========================================================

Interaction files are generated only when the corresponding interaction
is non-zero:

======================  ===========================================
File                    Interaction
======================  ===========================================
``coulombintra.def``    On-site Coulomb (intra-orbital)
``coulombinter.def``    Inter-site / inter-orbital Coulomb
``hund.def``            Hund coupling
``exchange.def``        Exchange interaction
``pairlift.def``        Pair-lift interaction
``pairhopp.def``        Pair-hopping interaction
``interall.def``        All remaining two-body interactions
======================  ===========================================

Green's function files are written when ``ioutputmode != 0``:

==================  =========================================================
File                Contents
==================  =========================================================
``greenone.def``    One-body Green's function measurement list
``greentwo.def``    Two-body Green's function measurement list
                    (not generated for UHF and H-wave)
==================  =========================================================

HPhi-Specific Files
^^^^^^^^^^^^^^^^^^^^

====================  =========================================================
File                  Contents
====================  =========================================================
``calcmod.def``       Calculation mode flags
``single.def``        Single-site observable list (non-pair output mode)
``pair.def``          Pair observable list (pair output mode)
``boost.def``         HPhi Boost-mode parameters (only when Boost is enabled)
``teone.def``         One-body time-evolution operators
                      (only when ``method = "timeevolution"``)
``tetwo.def``         Two-body time-evolution operators
                      (only when ``method = "timeevolution"``)
====================  =========================================================

mVMC-Specific Files
^^^^^^^^^^^^^^^^^^^^

======================  =========================================================
File                    Contents
======================  =========================================================
``gutzwilleridx.def``   Gutzwiller factor index table (one entry per site)
``jastrowidx.def``      Jastrow factor index table (one entry per site pair)
``orbitalidx.def``      Pairing-function (orbital) index table — used when
                        ``2Sz = 0`` or unspecified
``orbitalidxpara.def``  Spin-polarised orbital index table — used when
                        ``2Sz != 0`` or grand-canonical ensemble
``qptransidx.def``      Quantum-number projection / translation symmetry table
======================  =========================================================

Each index file maps variational parameters to lattice sites or site pairs and
records an optimisation flag (1 = optimised, 0 = fixed).

mVMC-Specific Input Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following keywords are written to ``modpara.def`` for mVMC and are ignored
by all other solvers.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Keyword
     - Description
   * - ``NVMCCalMode``
     - Calculation mode: ``0`` = variational optimisation,
       ``1`` = compute correlation functions with fixed wavefunction.
   * - ``NLanczosMode``
     - Power Lanczos correction: ``0`` = off, ``1`` = on.
   * - ``NDataIdxStart``
     - Start index of independent trial runs.
   * - ``NDataQtySmp``
     - Number of independent trial runs.
   * - ``NSPGaussLeg``
     - Number of Gauss-Legendre points for spin-quantum-number projection.
       Required when ``NSPStot`` is set.
   * - ``NSPStot``
     - Target total spin quantum number S for spin projection.
   * - ``NMPTrans``
     - Number of momentum/translation symmetry operations to project onto.
   * - ``NSROptItrStep``
     - Total number of stochastic reconfiguration (SR) optimisation steps.
   * - ``NSROptItrSmp``
     - Number of Monte Carlo samples per SR step.
   * - ``DSROptRedCut``
     - SR regularisation: eigenvalue cutoff (relative).
   * - ``DSROptStaDel``
     - SR regularisation: diagonal shift (stabiliser).
   * - ``DSROptStepDt``
     - SR optimisation step size (learning rate).
   * - ``NVMCWarmUp``
     - Number of Monte Carlo warm-up steps before sampling begins.
   * - ``NVMCInterval``
     - Interval (in MC steps) between successive samples.
   * - ``NVMCSample``
     - Total number of Monte Carlo samples per optimisation step.
   * - ``RndSeed``
     - Seed for the random number generator.
   * - ``NSplitSize``
     - Number of MPI groups for parallel trial splitting.
   * - ``NStore``
     - Flag to store intermediate wavefunction data to disk (``1`` = on).
   * - ``NSRCG``
     - Flag to use conjugate-gradient solver inside SR (``1`` = on).

Validation and Errors
---------------------

StdFace performs input validation during processing and reports errors to
``stdout``.

Common Error Messages
^^^^^^^^^^^^^^^^^^^^^

**ERROR ! <parameter> is NOT specified !**

A parameter expected by StdFace for the current model/lattice configuration
was not found in the input file. Add the indicated parameter to your input.

*Cause*: Integer parameters are initialized to a sentinel value (2147483647).
If unchanged after parsing, the parameter is considered missing.

**ERROR ! Keyword <parameter> is duplicated !**

The same keyword appears more than once in the input file.  StdFace does not
accept duplicate keys — it exits immediately with status ``-1``.

*Resolution*: Remove the duplicate entry, keeping only one occurrence.

**Check ! <parameter> is SPECIFIED but will NOT be USED.**

A parameter was specified in the input file that does not apply to the current
model/lattice combination.

*Resolution*: Remove or comment out the unused parameter from your input file.

*Cause*: Double parameters are initialized to NaN. If a non-NaN value is
detected for an inapplicable parameter, this error is raised.

What StdFace Validates
^^^^^^^^^^^^^^^^^^^^^^

Based on the error handling documentation:

- Presence of expected parameters (via sentinel value detection)
- Absence of inapplicable parameters (via NaN detection for doubles)
- Uniqueness of each keyword (duplicate keys are a fatal error)

What StdFace Does NOT Validate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following are unspecified or not validated by StdFace:

- Physical consistency of parameter values
- Numerical stability or correctness of the model
- Solver-specific parameter ranges or constraints
- Completeness of Wannier90 input files

Limitations and Non-Goals
-------------------------

StdFace is an input file generator, not a solver. The following are outside
its scope:

Solver-Specific Semantics
^^^^^^^^^^^^^^^^^^^^^^^^^

StdFace generates input files but does not interpret solver-specific options
beyond what is needed for file generation. Parameters like ``method`` and
``exct`` are passed through to output files; their validity is not checked
by StdFace.

Numerical Correctness
^^^^^^^^^^^^^^^^^^^^^

StdFace does not validate that the specified model parameters produce
physically meaningful or numerically stable results. Users are responsible
for ensuring appropriate parameter values.

Advanced Solver Options
^^^^^^^^^^^^^^^^^^^^^^^

Some solver options may not be exposed through StdFace's input format. For
advanced configuration, users may need to edit the generated files directly
or consult solver-specific documentation.

Output File Formats
^^^^^^^^^^^^^^^^^^^

The format and content of output files follow the input specifications of
each target solver.  For details of each file's internal format, refer to
the documentation of HPhi, mVMC, UHF, or H-wave respectively.
