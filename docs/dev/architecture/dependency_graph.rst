Dependency Graph
================

This document describes the dependencies between source files and components
in the StdFace codebase.

File Dependencies
-----------------

Header Inclusion Graph
^^^^^^^^^^^^^^^^^^^^^^

::

   dry.c
     └── version.h

   StdFace_main.c
     ├── StdFace_vals.h
     ├── StdFace_ModelUtil.h
     └── export_wannier90.h  [conditional: _HWAVE only]

   StdFace_ModelUtil.c
     ├── StdFace_vals.h
     ├── StdFace_ModelUtil.h
     └── setmemory.h

   ChainLattice.c
     ├── StdFace_vals.h
     └── StdFace_ModelUtil.h

   [All other lattice files follow the same pattern]

   setmemory.c
     └── setmemory.h

   export_wannier90.c
     ├── StdFace_vals.h
     └── export_wannier90.h

**Source reference**: Include statements in each ``.c`` file

Header Dependency Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^

::

   ┌─────────────────┐
   │   version.h     │  (standalone)
   └─────────────────┘

   ┌─────────────────┐
   │  setmemory.h    │  (standalone)
   └─────────────────┘

   ┌─────────────────┐
   │ StdFace_vals.h  │  (defines StdIntList)
   └────────┬────────┘
            │
            ▼
   ┌─────────────────┐     ┌─────────────────────┐
   │StdFace_ModelUtil.h│◄───│ export_wannier90.h  │
   │  (uses StdIntList) │    │ (conditional _HWAVE) │
   └─────────────────┘     └─────────────────────┘

**Key observations**:

- ``StdFace_vals.h`` is the foundational header (defines ``StdIntList``)
- ``StdFace_ModelUtil.h`` depends on ``StdFace_vals.h``
- ``export_wannier90.h`` depends on ``StdFace_vals.h`` (via forward declaration)
- ``version.h`` and ``setmemory.h`` are standalone

Call Graph
----------

Main Call Path
^^^^^^^^^^^^^^

::

   main()  [dry.c]
     │
     ├── printVersion()  [version.h]
     │
     └── StdFace_main()  [StdFace_main.c]
           │
           ├── StdFace_ResetVals()
           │
           ├── [Input parsing functions]
           │     ├── TrimSpaceQuote()
           │     ├── Text2Lower()
           │     └── StoreWithCheckDup_*()
           │
           ├── [Lattice constructor - one of:]
           │     ├── StdFace_Chain()        [ChainLattice.c]
           │     ├── StdFace_Ladder()       [Ladder.c]
           │     ├── StdFace_Tetragonal()   [SquareLattice.c]
           │     ├── StdFace_Triangular()   [TriangularLattice.c]
           │     ├── StdFace_Honeycomb()    [HoneycombLattice.c]
           │     ├── StdFace_Kagome()       [Kagome.c]
           │     ├── StdFace_Orthorhombic() [Orthorhombic.c]
           │     ├── StdFace_FCOrtho()      [FCOrtho.c]
           │     ├── StdFace_Pyrochlore()   [Pyrochlore.c]
           │     └── StdFace_Wannier90()    [Wannier90.c]
           │
           └── [Output functions - mode dependent]
                 ├── PrintLocSpin()
                 ├── PrintTrans()
                 ├── PrintInteractions()
                 ├── PrintModPara()
                 ├── PrintNamelist()
                 └── [mode-specific functions]

**Source reference**: ``src/dry.c``, ``src/StdFace_main.c``

Lattice Constructor Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each lattice constructor calls common utility functions from
``StdFace_ModelUtil.c``:

::

   StdFace_<Lattice>()
     │
     ├── StdFace_InitSite()       # Initialize site arrays
     │
     ├── StdFace_MallocInteractions()  # Allocate interaction arrays
     │
     ├── StdFace_SetLabel()       # Set site labels
     │
     ├── [Model-specific interaction setup]
     │     ├── StdFace_Hopping()      # Hopping terms
     │     ├── StdFace_HubbardLocal() # On-site U
     │     ├── StdFace_Coulomb()      # Inter-site V
     │     ├── StdFace_GeneralJ()     # Exchange J
     │     └── StdFace_MagField()     # Magnetic field
     │
     └── StdFace_PrintGeometry()  # Output geometry info

**Source reference**: Lattice constructor implementations

Memory Allocation Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Memory allocation functions in ``setmemory.c`` are called from:

::

   StdFace_MallocInteractions()  [StdFace_ModelUtil.c]
     │
     ├── i_1d_allocate()   # Integer arrays
     ├── i_2d_allocate()   # 2D integer arrays
     ├── d_1d_allocate()   # Double arrays
     ├── cd_1d_allocate()  # Complex double arrays
     └── cd_2d_allocate()  # 2D complex arrays

**Source reference**: ``src/StdFace_ModelUtil.c``, ``src/setmemory.c``

Component Responsibilities
--------------------------

+------------------------+------------------------------------------------+
| Component              | Responsibility                                 |
+========================+================================================+
| ``dry.c``              | CLI interface, argument handling               |
+------------------------+------------------------------------------------+
| ``StdFace_main.c``     | Input parsing, orchestration, output generation|
+------------------------+------------------------------------------------+
| ``StdFace_ModelUtil.c``| Lattice utilities, interaction setup           |
+------------------------+------------------------------------------------+
| ``<Lattice>.c`` files  | Lattice-specific geometry and neighbor logic   |
+------------------------+------------------------------------------------+
| ``setmemory.c``        | Memory allocation/deallocation                 |
+------------------------+------------------------------------------------+
| ``export_wannier90.c`` | Wannier90 format export (HWAVE only)           |
+------------------------+------------------------------------------------+

Build Dependencies
------------------

CMake Dependency Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   CMakeLists.txt (root)
     │
     ├── Sets project options (HPHI, MVMC, UHF, HWAVE)
     │
     └── add_subdirectory(src)
           │
           └── src/CMakeLists.txt
                 │
                 ├── Defines SOURCES_StdFace (common sources)
                 │
                 └── For each enabled mode:
                       ├── add_library(StdFace_<mode>)
                       ├── add_executable(<mode>_dry.out)
                       └── target_link_libraries(..., m)

**External dependencies**:

- Standard C library (``libc``)
- Math library (``libm``) - linked via ``-lm``

**Source reference**: ``CMakeLists.txt``, ``src/CMakeLists.txt``

Compile-Time Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^

The preprocessor defines create compile-time dependencies:

+-------------+-----------------------------------------------+
| Define      | Enables                                       |
+=============+===============================================+
| ``_HPhi``   | HPhi-specific keywords, output functions      |
+-------------+-----------------------------------------------+
| ``_mVMC``   | mVMC-specific keywords, variational outputs   |
+-------------+-----------------------------------------------+
| ``_UHF``    | UHF-specific keywords, mean-field outputs     |
+-------------+-----------------------------------------------+
| ``_HWAVE``  | HWAVE keywords, Wannier90 export functions    |
+-------------+-----------------------------------------------+

These are mutually exclusive in practice (each binary has one define).

**Source reference**: ``src/CMakeLists.txt:19-50``

Data Flow
---------

::

   Input File (stan.in)
         │
         ▼
   ┌─────────────────┐
   │ Parse keywords  │
   │ into StdIntList │
   └────────┬────────┘
            │
            ▼
   ┌─────────────────┐
   │ Lattice         │
   │ Constructor     │
   │                 │
   │ Populates:      │
   │ - trans[]       │
   │ - intr[]        │
   │ - transindx[][] │
   │ - intrindx[][]  │
   └────────┬────────┘
            │
            ▼
   ┌─────────────────┐
   │ Output          │
   │ Functions       │
   │                 │
   │ Write config    │
   │ files for       │
   │ target solver   │
   └────────┬────────┘
            │
            ▼
   Output Files (*.def, etc.)

The ``StdIntList`` structure serves as the central data container throughout
the entire processing pipeline.

**Source reference**: ``src/StdFace_vals.h`` (structure definition),
``src/StdFace_main.c`` (data flow implementation)
