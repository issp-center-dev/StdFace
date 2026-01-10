export_wannier90.h
==================

Purpose
-------

Wannier90 format export functions for HWAVE mode. These functions export geometry
and interaction data in a format compatible with Wannier90-style tight-binding
models.

.. important::

   This header is only available when ``_HWAVE`` is defined at compile time.
   The functions are wrapped in ``#if defined(_HWAVE)`` preprocessor conditionals.

**Source reference**: ``src/export_wannier90.h``

Public Types
------------

None.

Public Macros/Constants
-----------------------

.. c:macro:: EXPORT_WANNIER90_INCLUDED

   Include guard macro.

Public Functions
----------------

.. c:function:: void ExportGeometry(struct StdIntList *StdI)

   Exports lattice geometry information in Wannier90 format.

   :param StdI: Pointer to the configuration structure

   Writes geometry data including lattice vectors and site positions.
   Output format details are unspecified in the current code.

   **Availability**: Only when ``_HWAVE`` is defined.

.. c:function:: void ExportInteraction(struct StdIntList *StdI)

   Exports interaction parameters in Wannier90 format.

   :param StdI: Pointer to the configuration structure

   Writes hopping and interaction terms. Output format details are
   unspecified in the current code.

   **Availability**: Only when ``_HWAVE`` is defined.

Ownership and Lifetime Rules
----------------------------

- Functions do not allocate or free the ``StdI`` structure.
- Output files are created/written by these functions; file handles are
  managed internally.

Error Handling
--------------

Unspecified in the current code.

Thread / MPI Safety
-------------------

Unspecified in the current code.

Source Reference
----------------

- Header: ``src/export_wannier90.h``
- Implementation: ``src/export_wannier90.c``
- Build condition: ``_HWAVE`` defined (CMake option ``-DHWAVE=ON``)
