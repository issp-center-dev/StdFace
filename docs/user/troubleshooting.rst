Troubleshooting
===============

.. contents::
   :local:
   :depth: 2

Build Errors
------------

**No solver enabled**

If you run ``cmake ..`` without enabling any solver option, no executable will
be built.  Specify at least one solver:

.. code-block:: bash

   cmake .. -DHPHI=ON

**CMake version too old**

StdFace requires CMake 2.8.12 or later.  Check your version:

.. code-block:: bash

   cmake --version

Runtime Errors — Input File
----------------------------

**ERROR ! "=" is NOT found !**

A line in the input file does not contain ``=``.  Every non-comment line must
follow the ``key = value`` format.

*Resolution*: Check the indicated line for a missing ``=`` sign, or prefix it
with ``//`` to mark it as a comment.

**ERROR ! Keyword <x> is duplicated !**

The same keyword appears more than once in the input file.  StdFace does not
allow duplicate keys and exits immediately.

*Resolution*: Remove the duplicate entry, keeping only one occurrence.

**ERROR ! <x> is NOT specified !**

A parameter required for the current model/lattice combination is missing from
the input file.

*Resolution*: Add the indicated parameter.  Refer to :doc:`input_output` for
which parameters are required for each solver.

**Check ! <x> is SPECIFIED but will NOT be USED.**

A parameter was given that does not apply to the current model/lattice
combination.  StdFace continues but prints this warning.

*Resolution*: Remove or comment out the unused parameter to keep the input
file clean.

**ERROR ! Unsupported Keyword in Standard mode!**

The input file contains an unrecognised keyword.

*Resolution*: Check for typos.  Refer to :doc:`input_output` for the list of
supported keywords.

Runtime Errors — File Access
-----------------------------

**ERROR ! Cannot open input file <fname> !**

StdFace could not open the specified input file.

*Resolution*: Verify the file path and that the file exists:

.. code-block:: bash

   ls -la stan.in
   ./hphi_dry.out stan.in

When using ``lattice = "wannier90"``, also confirm that ``zvo_geom.dat``,
``zvo_hr.dat``, and ``zvo_ur.dat`` are present in the same directory.

Runtime Errors — Parameter Values
-----------------------------------

**ERROR ! Unsupported Solver : <method>**

The value given for the ``method`` keyword is not recognised by the current
solver.  This error is HPhi-specific.

*Resolution*: Use one of the supported values.  See :ref:`method-hphi` in
:doc:`input_output`.

**ERROR! (L, W, Height) and (a0W, ..., a2H) conflict !**

Both the shorthand dimension keywords (``W``, ``L``, ``Height``) and the
explicit lattice vector keywords (``a0W``, ``a1W``, …, ``a2H``) were
specified at the same time.

*Resolution*: Use one form or the other, not both.

**ERROR! <param> conflict !**

Two or more interaction parameters that cannot be used simultaneously were
both specified (for example, a scalar ``J`` and its tensor components ``Jx``,
``Jy``, ``Jz``).

*Resolution*: Use either the scalar shorthand or the individual components,
not a mixture.

MPI Note
--------

When StdFace is linked against an MPI library and exits normally, the MPI
runtime may print error messages such as::

   MPI_Finalize: MPI not initialized

StdFace anticipates this and prints the following line immediately before
calling ``exit()``:

.. code-block:: text

   #######  You DO NOT have to WORRY about the following MPI-ERROR MESSAGE.  #######

These MPI messages are harmless and can be safely ignored.
