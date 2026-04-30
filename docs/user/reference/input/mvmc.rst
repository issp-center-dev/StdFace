mVMC Parameters
===============

.. contents::
   :local:
   :depth: 2

Parameters specific to the mVMC solver (``mvmc_dry.out``).  For parameters
shared across solvers see :doc:`common`.

The following keywords are written to ``modpara.def`` and are ignored by all
other solvers.

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
