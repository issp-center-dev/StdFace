Directory Structure
===================

Repository Layout
-----------------

::

   StdFace/
   ├── CMakeLists.txt        # Top-level CMake configuration
   ├── README.md             # Project documentation
   ├── LICENSE               # GPL v3 license
   ├── src/                  # Source code AND headers
   │   ├── CMakeLists.txt    # Build rules for binaries/libraries
   │   ├── dry.c             # CLI entry point (main)
   │   ├── StdFace_main.c    # Core processing logic
   │   ├── StdFace_ModelUtil.c/.h   # Model utilities
   │   ├── StdFace_vals.h    # StdIntList parameter structure
   │   ├── setmemory.c/.h    # Memory allocation utilities
   │   ├── version.h         # Version macros
   │   ├── export_wannier90.c/.h  # Wannier90 export (HWAVE)
   │   ├── ChainLattice.c    # 1D chain lattice
   │   ├── Ladder.c          # 1D ladder lattice
   │   ├── SquareLattice.c   # 2D square/tetragonal lattice
   │   ├── TriangularLattice.c  # 2D triangular lattice
   │   ├── HoneycombLattice.c   # 2D honeycomb lattice
   │   ├── Kagome.c          # 2D Kagome lattice
   │   ├── Orthorhombic.c    # 3D orthorhombic lattice
   │   ├── FCOrtho.c         # 3D face-centered orthorhombic
   │   ├── Pyrochlore.c      # 3D pyrochlore lattice
   │   └── Wannier90.c       # Wannier90 interface
   ├── samples/              # Example input files
   │   └── hubbard/
   │       ├── default_model/
   │       │   └── stan.in
   │       └── wannier/
   │           ├── stan.in
   │           └── zvo_*.dat
   ├── test/                 # Test suite
   └── docs/                 # Documentation (Sphinx)

**Source reference**: Repository root directory listing

Header Files in src/
--------------------

.. important::

   StdFace does **not** follow the conventional C project layout of placing
   public headers in an ``include/`` directory. All headers reside in ``src/``
   alongside implementation files.

This design choice means:

- Headers and implementations are co-located
- No separate include path configuration needed beyond ``src/``
- The ``src/CMakeLists.txt`` includes ``include_directories(include)`` but this
  directory does not exist in the repository

**Source reference**: ``src/CMakeLists.txt:8``

Header Files
^^^^^^^^^^^^

+------------------------+--------------------------------------------------+
| Header                 | Purpose                                          |
+========================+==================================================+
| ``StdFace_ModelUtil.h``| Core model utilities (37 functions)              |
|                        |                                                  |
|                        | - Lattice constructors                           |
|                        | - Interaction generators                         |
|                        | - Input parsing helpers                          |
+------------------------+--------------------------------------------------+
| ``StdFace_vals.h``     | Parameter structure definition                   |
|                        |                                                  |
|                        | - ``struct StdIntList`` (~125 fields)            |
|                        | - Lattice geometry parameters                    |
|                        | - Model parameters (t, U, V, J, etc.)            |
+------------------------+--------------------------------------------------+
| ``setmemory.h``        | Memory allocation utilities (24 functions)       |
|                        |                                                  |
|                        | - Type-specific allocators (int, double, complex)|
|                        | - 1D, 2D, 3D array allocation                    |
+------------------------+--------------------------------------------------+
| ``export_wannier90.h`` | Wannier90 export (HWAVE mode only)               |
|                        |                                                  |
|                        | - ``ExportGeometry()``                           |
|                        | - ``ExportInteraction()``                        |
|                        | - Conditional: only when ``_HWAVE`` defined      |
+------------------------+--------------------------------------------------+
| ``version.h``          | Version information                              |
|                        |                                                  |
|                        | - ``printVersion()`` function                    |
|                        | - ``VERSION_MAJOR/MINOR/PATCH`` macros           |
+------------------------+--------------------------------------------------+

**Source reference**: ``docs/_meta/header_map.yaml``

Source Files
^^^^^^^^^^^^

**Core files** (compiled into all binary variants):

+------------------------+--------------------------------------------------+
| Source                 | Purpose                                          |
+========================+==================================================+
| ``dry.c``              | CLI entry point; argument parsing                |
+------------------------+--------------------------------------------------+
| ``StdFace_main.c``     | Main processing; input parsing; output generation|
+------------------------+--------------------------------------------------+
| ``StdFace_ModelUtil.c``| Implementation of model utility functions        |
+------------------------+--------------------------------------------------+
| ``setmemory.c``        | Memory allocation implementation                 |
+------------------------+--------------------------------------------------+
| ``export_wannier90.c`` | Wannier90 export implementation                  |
+------------------------+--------------------------------------------------+

**Lattice-specific files** (each implements one lattice type):

+------------------------+--------------------------------------------------+
| Source                 | Lattice Type                                     |
+========================+==================================================+
| ``ChainLattice.c``     | 1D chain                                         |
+------------------------+--------------------------------------------------+
| ``Ladder.c``           | 1D ladder (multi-leg)                            |
+------------------------+--------------------------------------------------+
| ``SquareLattice.c``    | 2D square/tetragonal                             |
+------------------------+--------------------------------------------------+
| ``TriangularLattice.c``| 2D triangular                                    |
+------------------------+--------------------------------------------------+
| ``HoneycombLattice.c`` | 2D honeycomb                                     |
+------------------------+--------------------------------------------------+
| ``Kagome.c``           | 2D Kagome                                        |
+------------------------+--------------------------------------------------+
| ``Orthorhombic.c``     | 3D simple orthorhombic/cubic                     |
+------------------------+--------------------------------------------------+
| ``FCOrtho.c``          | 3D face-centered orthorhombic/cubic              |
+------------------------+--------------------------------------------------+
| ``Pyrochlore.c``       | 3D pyrochlore                                    |
+------------------------+--------------------------------------------------+
| ``Wannier90.c``        | Wannier90 tight-binding import                   |
+------------------------+--------------------------------------------------+

**Source reference**: ``src/CMakeLists.txt:13-15``

samples/ Directory
------------------

Contains example input files demonstrating StdFace usage:

::

   samples/
   └── hubbard/
       ├── default_model/
       │   └── stan.in         # Basic Hubbard model on square lattice
       └── wannier/
           ├── stan.in         # Wannier90 mode example
           ├── zvo_geom.dat    # Geometry data
           ├── zvo_hr.dat      # Hopping data
           └── zvo_ur.dat      # Interaction data

The ``stan.in`` files use the standard ``keyword = value`` format:

.. code-block:: text

   model = "Hubbard"
   lattice = square
   W = 2
   L = 2
   t = 1
   U = 1

**Source reference**: ``samples/hubbard/default_model/stan.in``

test/ Directory
---------------

Contains test cases for validation. Test infrastructure is enabled by the
``TestStdFace`` CMake option (default: ``ON``).

**Source reference**: ``CMakeLists.txt:20-24``

Build Output
------------

After building with CMake, the following artifacts are produced (depending on
enabled options):

::

   build/
   ├── src/
   │   ├── hphi_dry.out       # HPhi mode binary (if -DHPHI=ON)
   │   ├── mvmc_dry.out       # mVMC mode binary (if -DMVMC=ON)
   │   ├── uhf_dry.out        # UHF mode binary (if -DUHF=ON)
   │   ├── hwave_dry.out      # HWAVE mode binary (if -DHWAVE=ON)
   │   ├── libStdFace_hphi.a  # HPhi static library
   │   ├── libStdFace_mvmc.a  # mVMC static library
   │   ├── libStdFace_uhf.a   # UHF static library
   │   └── libStdFace_hwave.a # HWAVE static library
   └── ...

**Source reference**: ``src/CMakeLists.txt:17-51``
