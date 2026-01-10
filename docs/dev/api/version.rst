version.h
=========

Purpose
-------

Version information and printing utility. Defines version macros following
Semantic Versioning and provides a function to print the version string.

**Source reference**: ``src/version.h``

Public Types
------------

None.

Public Macros/Constants
-----------------------

.. c:macro:: VERSION_MAJOR

   Major version number.

   **Current value**: 0

.. c:macro:: VERSION_MINOR

   Minor version number.

   **Current value**: 5

.. c:macro:: VERSION_PATCH

   Patch version number.

   **Current value**: 0

.. c:macro:: VERSION_PRERELEASE

   Prerelease identifier string (e.g., "alpha", "beta.1").

   **Current value**: "" (empty string)

.. c:macro:: _INCLUDE_VERSION

   Include guard macro.

Version Format
^^^^^^^^^^^^^^

The version follows `Semantic Versioning <http://semver.org>`_::

   <major>.<minor>.<patch>[-<prerelease>]

Examples:

- ``0.5.0`` (current release)
- ``1.0.0-alpha`` (hypothetical prerelease)

Public Functions
----------------

.. c:function:: void printVersion()

   Prints the version string to standard output.

   Output format::

      StdFace version <major>.<minor>.<patch>[-<prerelease>]

   The prerelease suffix is only printed if ``VERSION_PRERELEASE`` is non-empty.

   :returns: void

   **Implementation note**: This function is defined inline in the header file.

Ownership and Lifetime Rules
----------------------------

Not applicable (no dynamic allocation).

Error Handling
--------------

None (simple print function).

Thread / MPI Safety
-------------------

Uses ``printf()`` which is typically thread-safe for output, though interleaving
with other threads' output is possible. No MPI-specific considerations.

Source Reference
----------------

- Header: ``src/version.h``
- Called from: ``main()`` in ``src/dry.c``
