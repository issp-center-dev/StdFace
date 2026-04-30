H-wave Parameters
=================

.. contents::
   :local:
   :depth: 2

Parameters specific to the H-wave solver (``hwave_dry.out``).  For parameters
shared across all solvers see :doc:`common`.

Shared Parameters with UHF
---------------------------

The following parameters are common to both H-wave and UHF.  They are written
to ``modpara.def`` under the ``HWAVE_Cal_Parameters`` header.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Keyword
     - Default
     - Description
   * - ``Iteration_max``
     - 1000
     - Maximum number of self-consistent iterations.
   * - ``Mix``
     - 0.5
     - Mixing parameter for self-consistent update (0 < Mix ≤ 1).
   * - ``eps``
     - 8
     - Convergence threshold exponent; iteration stops when the residual
       is below 10\ :sup:`-eps`.
   * - ``EpsSlater``
     - 6
     - Convergence threshold exponent for the Slater determinant.
   * - ``RndSeed``
     - 123456789
     - Seed for the random number generator.
   * - ``NMPTrans``
     - 0
     - Number of momentum/translation symmetry operations.
   * - ``Wsub``, ``Lsub``, ``Hsub``
     - —
     - Sub-lattice dimensions (shorthand form).
   * - ``a0Wsub`` … ``a2Hsub``
     - —
     - Sub-lattice vectors (explicit form; cannot be combined with the
       shorthand form above).

H-wave-Specific Parameters
---------------------------

The following keywords are unique to H-wave and are not recognised by any
other solver.

.. list-table::
   :header-rows: 1
   :widths: 20 60 20

   * - Keyword
     - Description
     - Notes
   * - ``calcmode``
     - Calculation sub-mode.  ``"uhfr"`` selects real-space UHF;
       any other value (e.g. ``"uhfk"``, ``"rpa"``) selects k-space mode.
     - Required for H-wave.
   * - ``fileprefix``
     - Prefix string prepended to all output filenames in k-space mode
       (e.g. ``fileprefix = "run1"`` produces ``run1_geom.dat``).
     - k-space only.
   * - ``exportall``
     - When set to ``1``, zero-valued interaction elements are written to
       output files in k-space mode.  Default: omit zero elements.
     - k-space only.
   * - ``lattice_gp``
     - When set to ``0``, suppresses generation of ``lattice.gp``.
     -
