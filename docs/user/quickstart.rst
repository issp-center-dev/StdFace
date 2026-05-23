Quickstart
==========

This guide walks you through building StdFace and running your first input file
generation.

Prerequisites
-------------

- CMake 2.8.12 or later
- C99-compatible compiler (GCC, Clang, etc.)
- Standard math library (libm)

Building StdFace
----------------

StdFace uses CMake to build variant executables for different target solvers.
You must enable at least one solver mode during configuration.

1. Clone or download the StdFace source code.

2. Create a build directory and configure:

   .. code-block:: bash

      mkdir build && cd build
      cmake .. -DHPHI=ON

   This builds ``hphi_dry.out`` for the HPhi solver.

3. Compile:

   .. code-block:: bash

      make

4. The executable is created in the build directory.

Available Build Options
^^^^^^^^^^^^^^^^^^^^^^^

Select which solver(s) to build by enabling CMake options:

==============  ===================  =================================
Option          Executable           Target Solver
==============  ===================  =================================
``-DUHF=ON``    ``uhf_dry.out``      UHF (Unrestricted Hartree-Fock)
``-DMVMC=ON``   ``mvmc_dry.out``     mVMC (Variational Monte Carlo)
``-DHPHI=ON``   ``hphi_dry.out``     HPhi (Exact Diagonalization)
``-DHWAVE=ON``  ``hwave_dry.out``    H-wave
==============  ===================  =================================

You can enable multiple solvers in one build:

.. code-block:: bash

   cmake .. -DHPHI=ON -DMVMC=ON -DUHF=ON -DHWAVE=ON
   make

Running StdFace
---------------

Command Pattern
^^^^^^^^^^^^^^^

.. code-block:: bash

   ./<solver>_dry.out <input_file>

For example, with HPhi:

.. code-block:: bash

   ./hphi_dry.out stan.in

To display version information:

.. code-block:: bash

   ./hphi_dry.out -v

Input File
^^^^^^^^^^

StdFace reads a single input file containing model and lattice parameters in
``key = value`` format. The file name is arbitrary (commonly ``stan.in``).

Sample input file (from ``samples/hubbard/default_model/stan.in``):

.. code-block:: text

   model = "Hubbard"
   lattice = square
   W=2
   L=2
   method = "CG"
   2Sz = 0
   nelec = 4
   exct = 1
   t=1
   U=1

See :doc:`reference/index` for a full parameter reference.

Output Files
^^^^^^^^^^^^

StdFace generates input files for the target solver in the current working
directory.  The output files depend on the solver mode; see
:doc:`reference/output` for the complete list (e.g., ``namelist.def``,
``modpara.def``, ``trans.def``, and solver-specific files).

Your First Run
--------------

See :doc:`examples` for complete working examples including input files and
expected behaviour.  A minimal run looks like:

.. code-block:: bash

   ./hphi_dry.out stan.in

Next Steps
----------

- See :doc:`examples` for complete working examples
- See :doc:`reference/index` for detailed parameter reference
- See :doc:`troubleshooting` for build and runtime error resolution
