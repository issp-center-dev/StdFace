Wannier90 Import
================

Imports lattice geometry and interactions from external Wannier90 format
files instead of using a built-in lattice constructor.

==============  =====================================
Property        Value
==============  =====================================
Canonical name  ``wannier90``
Aliases         —
Dimensions      3 (geometry file always requires 3D lattice vectors; lower-
                dimensional systems can be represented with degenerate vectors)
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

Reader validation (C generators)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The geometry count selects the first ``NsiteUC`` orbitals. It must be positive
and no greater than ``num_wann`` in each matrix file that is read
(``_hr.dat``, ``_ur.dat``, ``_jr.dat``, and ``_dr.dat``). Additional orbitals
in a matrix file are validated but excluded from the generated model. Reducing
the geometry count may leave unused centre rows at the end of the geometry
file. A density file with fewer orbitals than the geometry is rejected rather
than silently filling the absent orbitals with zero.

Geometry vectors and the selected centres must each contain three finite real
values on one line. A matrix file must contain a comment line, positive
``num_wann`` and ``nrpts`` counts on separate lines, and exactly ``nrpts``
positive integer degeneracies (which may span multiple lines). The degeneracy
values retain the existing treatment; validation does not change their
physical weighting.

Each of the ``nrpts`` blocks contains exactly ``num_wann * num_wann`` records,
one per line, with seven fields. Every orbital pair in ``1..num_wann`` must
occur once in that block, all records must have the same R vector, and each R
block must be unique. The order of orbital pairs within a block is arbitrary.
R and -R are distinct blocks; the existing inversion-symmetry reduction still
applies.
Blank lines between or after records, CRLF line endings and a missing final
newline are accepted. Wrapped records, trailing text, extra fields or matrix
records, incomplete input and malformed integer tokens are errors. Counts and
allocation products must fit the implementation's integer and memory limits;
matrix record counts are limited to ``INT_MAX``. R components must fit ``int``
and cannot be ``INT_MIN``.

The real-space cutoffs (``cutoff_length_t``, ``cutoff_length_U`` and
``cutoff_length_J``) use the centres of the orbitals named in columns 4 and 5
of each record: the orbital in the home cell and the orbital in cell R.
The result therefore does not depend on the record order. Earlier C
generators used the position of the record in its block instead. That gives
the same result for files written with column 5 varying fastest, such as
RESPACK output. For Wannier90's own ``seedname_hr.dat``, which varies column 4
fastest, the generated terms can change when a positive cutoff length is used,
the two orbitals have different centres, and R is nonzero: the old and corrected
distances may fall on opposite sides of the cutoff. With distance cutoffs
disabled, this correction does not change the generated terms.

Reader errors identify the input file and read position on ``stderr`` and
stop generation. These structural checks do not establish physical validity
or change the handling of nonfinite matrix values. The Python reader is
unchanged: it retains the old distance calculation and can disagree with the
corrected C generators for Wannier90-order inputs with distance cutoffs.
Selecting fewer orbitals through the geometry count can still raise
``IndexError`` in the Python reader.

Density R domain
^^^^^^^^^^^^^^^^

With ``doublecounting`` enabled, the density file need not list every R vector
of the interaction files. The generators take the density at an R vector
that ``_dr.dat`` does not list as zero, whether or not that R lies between the
smallest and largest R of the file, just as an R vector that ``_hr.dat``,
``_ur.dat`` or ``_jr.dat`` does not list has no term. The Hartree corrections
use only the density at R = 0 and are not affected. For an unlisted R,
``initial.def`` gets no off-diagonal element, and ``doublecounting = "full"``
adds no exchange (Fock) correction for the U and J terms at that R. No warning
is printed, so list every R whose density should enter the correction.

Earlier C generators read outside the density array when an interaction term
had an R beyond the smallest or largest R of the density file, and the Python
reader raised ``KeyError`` for such an R. Both now follow the rule above.

File Relationships
^^^^^^^^^^^^^^^^^^

.. code-block:: text

   stan.in  (lattice = "wannier90")
       |
       +-- zvo_geom.dat  (lattice geometry)
       +-- zvo_hr.dat    (hopping t)
       +-- zvo_ur.dat    (interaction U)
