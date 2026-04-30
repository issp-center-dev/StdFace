Triangular Lattice
==================

Two-dimensional triangular lattice.

.. figure:: /_static/lattice/triangular.png
   :width: 500px
   :align: center

   Triangular lattice.  Three nearest-neighbor bond directions:
   along W (magenta/red, t\ :sub:`0`, V\ :sub:`0`, J\ :sub:`0`),
   along L (green, t\ :sub:`1`, V\ :sub:`1`, J\ :sub:`1`),
   and along W−L (blue, t\ :sub:`2`, V\ :sub:`2`, J\ :sub:`2`).
   Primed labels (dashed) indicate 2nd nearest-neighbor bonds;
   double-primed (dotted) indicate 3rd nearest-neighbor bonds.

==============  =====================================
Property        Value
==============  =====================================
Canonical name  ``triangular``
Aliases         ``triangularlattice``
Dimensions      2D
Boost mode      Not supported
==============  =====================================

Interactions
------------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Bond class
     - Directions
     - Keywords
   * - Nearest neighbor — along W
     - 0
     - ``t0``, ``V0``, ``J0``
   * - Nearest neighbor — along L
     - 1
     - ``t1``, ``V1``, ``J1``
   * - Nearest neighbor — along W−L
     - 2
     - ``t2``, ``V2``, ``J2``
   * - 2nd nearest neighbor
     - 0, 1, 2
     - ``t0'``, ``V0'``, ``J0'``; ``t1'``, ``V1'``, ``J1'``; ``t2'``, ``V2'``, ``J2'``
   * - 3rd nearest neighbor
     - 0, 1, 2
     - ``t0''``, ``V0''``, ``J0''``; ``t1''``, ``V1''``, ``J1''``; ``t2''``, ``V2''``, ``J2''``
   * - On-site
     - —
     - ``U``, ``mu``, ``h``, ``Gamma``, ``Gamma_y``, ``D``

The isotropic shorthand ``t`` / ``V`` / ``J`` sets all three nearest-neighbor
directions simultaneously.
