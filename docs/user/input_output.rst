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

The specific output files depend on the solver mode (determined at compile time)
and are unspecified in the current code documentation.

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

Comment syntax is unspecified in the current code documentation. The sample
files do not contain comments.

Ordering
^^^^^^^^

Key ordering within the input file is unspecified. The sample files do not
suggest any required ordering.

Common Keys (Observed in Samples)
---------------------------------

The following keys are observed in ``samples/hubbard/default_model/stan.in``:

===========  =============  ================================================
Key          Sample Value   Description
===========  =============  ================================================
``model``    ``"Hubbard"``  Physical model type (observed in sample)
``lattice``  ``square``     Lattice geometry (observed in sample)
``W``        ``2``          Lattice dimension in W direction (observed)
``L``        ``2``          Lattice dimension in L direction (observed)
``method``   ``"CG"``       Solver method (solver-specific, observed)
``2Sz``      ``0``          Total spin projection 2*Sz (observed in sample)
``nelec``    ``4``          Number of electrons (observed in sample)
``exct``     ``1``          Number of excited states (solver-specific)
``t``        ``1``          Hopping parameter (observed in sample)
``U``        ``1``          On-site Coulomb interaction (observed in sample)
===========  =============  ================================================

**Evidence**: ``samples/hubbard/default_model/stan.in``

Key Details
^^^^^^^^^^^

``model``
   Physical model type. Observed value: ``"Hubbard"``. Other model types
   (e.g., spin, Kondo) are referenced in dev docs but not observed in samples.

``lattice``
   Lattice geometry. Observed values:

   - ``square`` - 2D square lattice (default_model sample)
   - ``"wannier90"`` - Import from Wannier90 files (wannier sample)

   Other lattice types documented in dev docs: ``chain``, ``ladder``,
   ``triangular``, ``honeycomb``, ``kagome``, ``tetragonal``, ``orthorhombic``,
   ``fcortho``, ``pyrochlore``.

``W``, ``L``
   Lattice dimensions. Interpretation depends on lattice type.

``method``
   Solver method. Observed value: ``"CG"`` (Conjugate Gradient). This key
   appears solver-specific; its valid values depend on the target solver.

``2Sz``
   Total spin projection (2 times Sz). Used in quantum number specification.

``nelec``
   Number of electrons in the system.

``exct``
   Number of excited states to compute. Appears solver-specific (e.g., for
   exact diagonalization).

``t``
   Nearest-neighbor hopping amplitude. Present in default_model sample but
   absent in wannier sample (hopping comes from external files).

``U``
   On-site Coulomb interaction strength. Present in default_model sample but
   absent in wannier sample (interaction comes from external files).

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

StdFace generates input files for the target solver in the current working
directory.

**Important**: The specific output file names, formats, and contents are
solver-dependent and unspecified in the current code documentation. Users
should consult the generated files after running StdFace.

The solver mode is determined at compile time:

==============  ===================  =====================================
Build Option    Executable           Target Solver
==============  ===================  =====================================
``-DUHF=ON``    ``uhf_dry.out``      UHF (Unrestricted Hartree-Fock)
``-DMVMC=ON``   ``mvmc_dry.out``     mVMC (Variational Monte Carlo)
``-DHPHI=ON``   ``hphi_dry.out``     HPhi (Exact Diagonalization)
``-DHWAVE=ON``  ``hwave_dry.out``    H-wave
==============  ===================  =====================================

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

The format and content of output files are determined by the solver mode and
are unspecified in StdFace's documentation. Users should refer to the
documentation of the target solver (HPhi, mVMC, UHF, or H-wave) for details.
