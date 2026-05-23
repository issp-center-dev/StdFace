Build System and Options
========================

This document describes the CMake build system and compile-time options for
StdFace.

Overview
--------

StdFace uses CMake as its build system. The build produces variant executables
for different target solvers (HPhi, mVMC, UHF, HWAVE), each with mode-specific
functionality enabled via preprocessor defines.

**Source references**:

- ``CMakeLists.txt`` (root)
- ``src/CMakeLists.txt``

Build Requirements
------------------

- CMake 2.8.12 or later
- C99-compatible compiler
- Standard math library (``libm``)
- Optional: MPI library (for MPI-aware error handling)

CMake Options
-------------

The following CMake options control which executables are built:

============  ==============  =======================================
Option        Default         Description
============  ==============  =======================================
``UHF``       OFF             Build ``uhf_dry.out`` for UHF solver
``MVMC``      OFF             Build ``mvmc_dry.out`` for mVMC solver
``HPHI``      OFF             Build ``hphi_dry.out`` for HPhi solver
``HWAVE``     OFF             Build ``hwave_dry.out`` for HWAVE solver
============  ==============  =======================================

**Source reference**: ``CMakeLists.txt`` lines 8-11

.. code-block:: cmake

   option(UHF "Build uhf_dry.out" OFF)
   option(MVMC "Build mvmc_dry.out" OFF)
   option(HPHI "Build hphi_dry.out" OFF)
   option(HWAVE "Build hwave_dry.out" OFF)

Build Configuration Examples
----------------------------

Build for HPhi solver:

.. code-block:: bash

   mkdir build && cd build
   cmake .. -DHPHI=ON
   make

Build for multiple solvers:

.. code-block:: bash

   cmake .. -DHPHI=ON -DMVMC=ON -DUHF=ON -DHWAVE=ON
   make

Preprocessor Defines
--------------------

Each build mode defines a preprocessor macro that enables mode-specific code:

=============  ====================  ==========================================
CMake Option   Preprocessor Define   Affected Features
=============  ====================  ==========================================
``UHF``        ``_UHF``              UHF-specific fields in StdIntList
``MVMC``       ``_mVMC``             mVMC orbital/variational functions
``HPHI``       ``_HPhi``             HPhi Lanczos and boost functions
``HWAVE``      ``_HWAVE``            HWAVE export functions (Wannier90 export)
=============  ====================  ==========================================

**Source reference**: ``src/CMakeLists.txt``

.. code-block:: cmake

   if(UHF)
      add_library(StdFace_uhf STATIC ${SOURCES_StdFace})
      target_compile_definitions(StdFace_uhf PUBLIC _UHF)
      add_executable(uhf_dry.out dry.c)
      target_compile_definitions(uhf_dry.out PUBLIC _UHF)
      target_link_libraries(uhf_dry.out PUBLIC StdFace_uhf m)
   endif()

Build Artifacts
---------------

Each enabled mode produces:

1. A static library: ``StdFace_<mode>``
2. An executable: ``<mode>_dry.out``

===========  ====================  ==================
Mode         Library               Executable
===========  ====================  ==================
UHF          ``StdFace_uhf``       ``uhf_dry.out``
MVMC         ``StdFace_mvmc``      ``mvmc_dry.out``
HPHI         ``StdFace_hphi``      ``hphi_dry.out``
HWAVE        ``StdFace_hwave``     ``hwave_dry.out``
===========  ====================  ==================

Conditional Compilation
-----------------------

Mode-Specific Code Blocks
^^^^^^^^^^^^^^^^^^^^^^^^^

Mode-specific code is wrapped in preprocessor conditionals:

.. code-block:: c

   #if defined(_HPhi)
   void StdFace_Chain_Boost(struct StdIntList *StdI);
   void StdFace_Ladder_Boost(struct StdIntList *StdI);
   void StdFace_Honeycomb_Boost(struct StdIntList *StdI);
   void StdFace_Kagome_Boost(struct StdIntList *StdI);
   #endif

   #if defined(_mVMC)
   void StdFace_generate_orb(struct StdIntList *StdI);
   void StdFace_Proj(struct StdIntList *StdI);
   void PrintJastrow(struct StdIntList *StdI);
   #endif

**Source reference**: ``src/StdFace_ModelUtil.h``

export_wannier90.h (HWAVE Only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``export_wannier90.h`` header provides Wannier90 export functions that are
only available when ``_HWAVE`` is defined:

.. code-block:: c

   #if defined(_HWAVE)
   void ExportGeometry(struct StdIntList *StdI);
   void ExportInteraction(struct StdIntList *StdI);
   #endif

**Source reference**: ``src/export_wannier90.h``

These functions export lattice geometry and interaction data in Wannier90
format for downstream processing by the HWAVE solver.

StdIntList Mode-Specific Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``StdIntList`` structure contains additional fields depending on build mode:

- **HPhi mode** (``_HPhi``): Lanczos parameters, spectrum calculation settings,
  time evolution parameters
- **mVMC mode** (``_mVMC``): Variational parameters, sublattice settings,
  orbital indices
- **UHF mode** (``_UHF``): Mean-field iteration parameters, sublattice settings
- **HWAVE mode** (``_HWAVE``): Wannier90 export settings (similar to UHF)

**Source reference**: ``src/StdFace_vals.h``

Source Files
------------

The following source files are compiled into each StdFace variant:

==========================  =============================================
File                        Description
==========================  =============================================
``StdFace_main.c``          Main entry and orchestration
``StdFace_ModelUtil.c``     Core model utilities
``ChainLattice.c``          1D chain lattice
``Ladder.c``                Ladder lattice
``FCOrtho.c``               Face-centered orthorhombic lattice
``HoneycombLattice.c``      Honeycomb lattice
``Kagome.c``                Kagome lattice
``Orthorhombic.c``          Orthorhombic/cubic lattice
``Pyrochlore.c``            Pyrochlore lattice
``SquareLattice.c``         Square lattice
``TriangularLattice.c``     Triangular lattice
``Wannier90.c``             Wannier90 import
``export_wannier90.c``      Wannier90 export (HWAVE)
``version.c``               Version information
``setmemory.c``             Memory allocation utilities
==========================  =============================================

**Source reference**: ``src/CMakeLists.txt``

Shared Entry Point
------------------

All variants share the same entry point code in ``dry.c``:

.. code-block:: c

   int main(int argc, char* argv[]) {
       StdFace_main(argv[1]);
       return 0;
   }

The mode-specific behavior is determined by the preprocessor defines, not by
runtime arguments.

Source References
-----------------

- CMake configuration: ``CMakeLists.txt``, ``src/CMakeLists.txt``
- Conditional declarations: ``src/StdFace_ModelUtil.h``
- Mode-specific structure fields: ``src/StdFace_vals.h``
- HWAVE export functions: ``src/export_wannier90.h``
