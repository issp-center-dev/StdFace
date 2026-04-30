UHF Parameters
==============

Parameters specific to the UHF solver (``uhf_dry.out``).  For parameters
shared across all solvers see :doc:`common`.

The following parameters are common to both UHF and H-wave.  They are written
to ``modpara.def`` under the ``UHF_Cal_Parameters`` header.

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
