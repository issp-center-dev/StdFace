Output Files
============

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

.. contents::
   :local:
   :depth: 2

Files Common to All Solvers
----------------------------

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
--------------------

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
--------------------

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
