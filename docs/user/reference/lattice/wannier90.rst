Wannier90 Import
================

Imports lattice geometry and interactions from external Wannier90 format
files instead of using a built-in lattice constructor.

==============  =====================================
Property        Value
==============  =====================================
Canonical name  ``wannier90``
Aliases         —
Dimensions      Any
Boost mode      Not supported
==============  =====================================

Required Input Files
--------------------

When ``lattice = "wannier90"`` is specified, three additional files must be
present in the same directory as ``stan.in``:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - File
     - Contents
   * - ``zvo_geom.dat``
     - Lattice vectors and atomic positions
   * - ``zvo_hr.dat``
     - Hopping (transfer) integrals in Wannier90 format
   * - ``zvo_ur.dat``
     - On-site interaction parameters in Wannier90 format

Example ``stan.in``
-------------------

.. code-block:: text

   model   = "Hubbard"
   lattice = "wannier90"
   W       = 2
   L       = 2
   method  = "CG"
   2Sz     = 0
   nelec   = 4
   exct    = 1

Note: ``t`` and ``U`` are absent; hopping and interaction values come from
the Wannier90 files.

File Formats
------------

Geometry File (``zvo_geom.dat``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

     1.0000000000   0.0000000000   0.0000000000
     0.0000000000   1.0000000000   0.0000000000
     0.0000000000   0.0000000000   1.0000000000
     1
     0.5000000000 0.5000000000 0.5000000000

- Lines 1–3: Lattice vectors (3×3 matrix, one vector per line)
- Line 4: Number of atoms in the unit cell
- Lines 5+: Fractional coordinates of each atom

Hopping File (``zvo_hr.dat``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   wannier90 format for vmcdry.out or HPhi -sdry
            1
            9
       1    1    1    1    1    1    1    1    1
      -1   -1    0    1    1    0.0  0.0
      -1    0    0    1    1   -1.0  0.0
      ...

- Line 1: Header comment
- Line 2: Number of Wannier functions
- Line 3: Number of R vectors
- Line 4: Degeneracy weights for each R vector
- Lines 5+: ``R_x  R_y  R_z  orbital_i  orbital_j  Re(t)  Im(t)``

Interaction File (``zvo_ur.dat``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Same format as ``zvo_hr.dat``.  The entry at R = (0, 0, 0) with value 1.0
represents the on-site Hubbard U interaction.

File Relationships
^^^^^^^^^^^^^^^^^^

.. code-block:: text

   stan.in  (lattice = "wannier90")
       |
       +-- zvo_geom.dat  (lattice geometry)
       +-- zvo_hr.dat    (hopping t)
       +-- zvo_ur.dat    (interaction U)
