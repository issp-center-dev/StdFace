Overview
========

Purpose
-------

StdFace is an input file generator for quantum lattice model solvers. It reads a
simple, human-readable input file (``stan.in`` format) and generates the detailed
configuration files required by target solvers:

- **HPhi**: Exact diagonalization solver
- **mVMC**: Many-variable Variational Monte Carlo solver
- **UHF**: Unrestricted Hartree-Fock solver
- **H-wave**: H-wave quantum lattice solver

The "dry run" naming convention (e.g., ``hphi_dry.out``) indicates that StdFace
generates input files without executing the solver itself.

**Source reference**: ``README.md``, ``src/StdFace_main.c:2-33``

Design Philosophy
-----------------

StdFace follows a **single-source, multi-target** design:

1. **Unified Input Format**: Users specify lattice geometry and model parameters
   in a simple ``keyword = value`` format.

2. **Build-Time Mode Selection**: The target solver is determined at compile time
   via preprocessor defines (``_HPhi``, ``_mVMC``, ``_UHF``, ``_HWAVE``). This
   produces four separate binaries from the same source code.

3. **Lattice Abstraction**: Lattice-specific logic is encapsulated in dedicated
   constructor functions (e.g., ``StdFace_Chain``, ``StdFace_Honeycomb``), enabling
   support for multiple lattice geometries.

**Source reference**: ``src/CMakeLists.txt:17-51``, ``src/StdFace_main.c:2915-2945``

Header Location Note
--------------------

.. important::

   Unlike typical C projects, StdFace does **not** have an ``include/`` directory.
   All header files are located in ``src/``:

   - ``src/StdFace_ModelUtil.h`` - Core model utilities (37 functions)
   - ``src/StdFace_vals.h`` - Parameter structure ``StdIntList`` (~125 fields)
   - ``src/setmemory.h`` - Memory allocation utilities (24 functions)
   - ``src/export_wannier90.h`` - Wannier90 export (HWAVE-only, 2 functions)
   - ``src/version.h`` - Version information (1 function, 4 macros)

   See :doc:`directory_structure` for details.

Core Components
---------------

Entry Point
^^^^^^^^^^^

The CLI entry point is ``main()`` in ``src/dry.c``, which:

1. Parses command-line arguments
2. Handles ``-v`` for version display
3. Dispatches to ``StdFace_main(fname)`` for processing

**Source reference**: ``src/dry.c:34-47``

Main Processing
^^^^^^^^^^^^^^^

``StdFace_main()`` in ``src/StdFace_main.c`` performs:

1. **Initialization**: Allocates and resets ``StdIntList`` structure
2. **Input Parsing**: Reads ``keyword = value`` pairs from input file
3. **Lattice Construction**: Calls appropriate lattice function based on input
4. **Output Generation**: Writes solver-specific configuration files

**Source reference**: ``src/StdFace_main.c:2454-3066``

Supported Lattices
^^^^^^^^^^^^^^^^^^

StdFace supports these lattice types (each has a dedicated constructor function):

+------------------+-------------------------+----------------------------+
| Dimension        | Lattice                 | Function                   |
+==================+=========================+============================+
| 1D               | Chain                   | ``StdFace_Chain``          |
+------------------+-------------------------+----------------------------+
| 1D               | Ladder                  | ``StdFace_Ladder``         |
+------------------+-------------------------+----------------------------+
| 2D               | Square/Tetragonal       | ``StdFace_Tetragonal``     |
+------------------+-------------------------+----------------------------+
| 2D               | Triangular              | ``StdFace_Triangular``     |
+------------------+-------------------------+----------------------------+
| 2D               | Honeycomb               | ``StdFace_Honeycomb``      |
+------------------+-------------------------+----------------------------+
| 2D               | Kagome                  | ``StdFace_Kagome``         |
+------------------+-------------------------+----------------------------+
| 3D               | Orthorhombic            | ``StdFace_Orthorhombic``   |
+------------------+-------------------------+----------------------------+
| 3D               | Face-Centered Ortho.    | ``StdFace_FCOrtho``        |
+------------------+-------------------------+----------------------------+
| 3D               | Pyrochlore              | ``StdFace_Pyrochlore``     |
+------------------+-------------------------+----------------------------+
| N/A              | Wannier90 import        | ``StdFace_Wannier90``      |
+------------------+-------------------------+----------------------------+

**Source reference**: ``src/StdFace_main.c:23-33``, ``src/StdFace_ModelUtil.h``

Build-Time Mode Selection
-------------------------

The same source code produces different binaries depending on CMake options:

+----------------+--------------+-------------------+------------------------+
| CMake Option   | Define       | Binary            | Library                |
+================+==============+===================+========================+
| ``-DHPHI=ON``  | ``_HPhi``    | ``hphi_dry.out``  | ``libStdFace_hphi.a``  |
+----------------+--------------+-------------------+------------------------+
| ``-DMVMC=ON``  | ``_mVMC``    | ``mvmc_dry.out``  | ``libStdFace_mvmc.a``  |
+----------------+--------------+-------------------+------------------------+
| ``-DUHF=ON``   | ``_UHF``     | ``uhf_dry.out``   | ``libStdFace_uhf.a``   |
+----------------+--------------+-------------------+------------------------+
| ``-DHWAVE=ON`` | ``_HWAVE``   | ``hwave_dry.out`` | ``libStdFace_hwave.a`` |
+----------------+--------------+-------------------+------------------------+

All options default to ``OFF``. One or more may be enabled simultaneously.

**Source reference**: ``CMakeLists.txt:4-14``, ``src/CMakeLists.txt:17-51``

Thread and MPI Safety
---------------------

Unspecified in the current code. The codebase does not appear to use explicit
threading or MPI constructs in the ``src/`` files examined.
