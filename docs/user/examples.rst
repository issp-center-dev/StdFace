Examples
========

This section provides complete working examples for different use cases.

Hubbard Model (Default)
-----------------------

This example demonstrates a basic Hubbard model on a 2x2 square lattice using
the standard lattice definition. This is the recommended starting point.

**Source**: ``samples/hubbard/default_model/``

Input File
^^^^^^^^^^

Create ``stan.in``:

.. code-block:: text

   model = "Hubbard"
   lattice = square
   W = 2
   L = 2
   method = "CG"
   2Sz = 0
   nelec = 4
   exct = 1
   t = 1
   U = 1

See :doc:`input_output` for a full description of each parameter.

Running the Example
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cd samples/hubbard/default_model
   ../../../build/hphi_dry.out stan.in

Replace ``hphi_dry.out`` with the appropriate executable for your target solver.

Expected Behavior
^^^^^^^^^^^^^^^^^

StdFace reads the input file, constructs a 2x2 square lattice Hubbard model,
and generates solver input files in the current directory. The specific output
files depend on the solver mode used.

Hubbard Model (Wannier90-based)
-------------------------------

This advanced example demonstrates importing a tight-binding model from
Wannier90 format files. This approach is useful when you have an existing
Wannier90 tight-binding Hamiltonian.

**Source**: ``samples/hubbard/wannier/``

Required Files
^^^^^^^^^^^^^^

This example requires multiple input files:

1. ``stan.in`` - Main StdFace input file
2. ``zvo_geom.dat`` - Geometry data (lattice vectors and atomic positions)
3. ``zvo_hr.dat`` - Hopping parameters in Wannier90 format
4. ``zvo_ur.dat`` - On-site interaction parameters in Wannier90 format

Main Input File (stan.in)
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   model = "Hubbard"
   lattice = "wannier90"
   W = 2
   L = 2
   method = "CG"
   2Sz = 0
   nelec = 4
   exct = 1

Key difference from the default example: ``lattice = "wannier90"`` instructs
StdFace to read lattice and hopping information from Wannier90 format files
instead of using built-in lattice constructors.

Geometry File (zvo_geom.dat)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file defines the lattice vectors and atomic positions:

.. code-block:: text

     1.0000000000   0.0000000000   0.0000000000
     0.0000000000   1.0000000000   0.0000000000
     0.0000000000   0.0000000000   1.0000000000
     1
     0.5000000000 0.5000000000 0.5000000000

Format:

- Lines 1-3: Lattice vectors (3x3 matrix, one vector per line)
- Line 4: Number of atoms in the unit cell
- Line 5+: Fractional coordinates of each atom

Hopping File (zvo_hr.dat)
^^^^^^^^^^^^^^^^^^^^^^^^^

This file defines the hopping (transfer) integrals in Wannier90 format:

.. code-block:: text

   wannier90 format for vmcdry.out or HPhi -sdry
            1
            9
       1    1    1    1    1    1    1    1    1
      -1   -1    0    1    1    0.0  0.0
      -1    0    0    1    1   -1.0  0.0
      ...

Format:

- Line 1: Header comment
- Line 2: Number of Wannier functions
- Line 3: Number of R vectors
- Line 4: Degeneracy weights for each R vector
- Lines 5+: Hopping data (R_x, R_y, R_z, orbital_i, orbital_j, Re(t), Im(t))

Interaction File (zvo_ur.dat)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file defines the on-site interaction in the same format as the hopping
file:

.. code-block:: text

   wannier90 format for vmcdry.out or HPhi -sdry
            1
            9
       1    1    1    1    1    1    1    1    1
      -1   -1    0    1    1    0.0  0.0
      ...
       0    0    0    1    1    1.0  0.0
      ...

The entry at R = (0, 0, 0) with value 1.0 represents the on-site Hubbard U.

File Relationships
^^^^^^^^^^^^^^^^^^

.. code-block:: text

   stan.in
     |
     +-- lattice = "wannier90"
           |
           +-- zvo_geom.dat  (geometry)
           +-- zvo_hr.dat    (hopping t)
           +-- zvo_ur.dat    (interaction U)

When ``lattice = "wannier90"`` is specified, StdFace automatically looks for
these Wannier90 format files in the same directory as the input file.

Running the Example
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cd samples/hubbard/wannier
   ../../../build/hphi_dry.out stan.in

Ensure all four files (``stan.in``, ``zvo_geom.dat``, ``zvo_hr.dat``,
``zvo_ur.dat``) are present in the working directory.

When to Use Wannier90 Format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use Wannier90-based input when:

- You have an existing tight-binding model from DFT+Wannier90 calculations
- You need custom hopping patterns not supported by built-in lattices
- You want to import realistic band structures from first-principles calculations

For simple models with standard lattice geometries, the default approach
(built-in lattice constructors) is recommended.

Lattice Types Reference
-----------------------

See :ref:`Lattice Types <lattice-types>` in :doc:`input_output` for a full
list of supported lattice geometries.
