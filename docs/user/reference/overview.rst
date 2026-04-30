Overview
========

StdFace is an input file generator for quantum lattice model solvers. It reads
a configuration file specifying the physical model and lattice geometry, then
generates solver-specific input files.

- **Input**: StdFace reads a single text configuration file (commonly named
  ``stan.in``) containing model and lattice parameters.
- **Output**: StdFace writes one or more input files for the target solver
  (HPhi, mVMC, UHF, or H-wave) in the current working directory.

The specific output files depend on the solver mode, which is determined at
compile time.  See :doc:`output` for a complete list.
