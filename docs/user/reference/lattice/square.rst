Square Lattice
==============

Two-dimensional square (tetragonal) lattice.

.. figure:: /_static/lattice/square.png
   :width: 380px
   :align: center

   Square lattice.  Solid (magenta/red): nearest-neighbor bonds
   along W (t\ :sub:`0`, V\ :sub:`0`, J\ :sub:`0`) and
   along L (t\ :sub:`1`, V\ :sub:`1`, J\ :sub:`1`).
   Dashed (orange/cyan): 2nd nearest-neighbor diagonal bonds
   (t\ :sub:`0`\ ', t\ :sub:`1`\ ').
   Dotted (blue): 3rd nearest-neighbor bonds
   (t\ :sub:`0`\ '', t\ :sub:`1`\ '').

==============  =====================================
Property        Value
==============  =====================================
Canonical name  ``tetragonal``
Aliases         ``square``, ``squarelattice``
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
   * - 2nd nearest neighbor — diagonal (W+L)
     - 0
     - ``t0'``, ``V0'``, ``J0'``
   * - 2nd nearest neighbor — anti-diagonal (W−L)
     - 1
     - ``t1'``, ``V1'``, ``J1'``
   * - 3rd nearest neighbor — along 2W
     - 0
     - ``t0''``, ``V0''``, ``J0''``
   * - 3rd nearest neighbor — along 2L
     - 1
     - ``t1''``, ``V1''``, ``J1''``
   * - On-site
     - —
     - ``U``, ``mu``, ``h``, ``Gamma``, ``Gamma_y``, ``D``

When W = L (isotropic), the shorthands ``t`` / ``V`` / ``J``
set both directions simultaneously (equivalent to ``t0`` = ``t1``).
