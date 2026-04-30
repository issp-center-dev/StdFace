Ladder Lattice
==============

One-dimensional two-leg (or multi-leg) ladder.

.. figure:: /_static/lattice/ladder.png
   :width: 420px
   :align: center

   Ladder lattice.  Solid (magenta): rung bonds
   (t\ :sub:`0`, V\ :sub:`0`, J\ :sub:`0`).
   Dashed (green): leg bonds
   (t\ :sub:`2`, V\ :sub:`2`, J\ :sub:`2`).
   Dashed (blue/cyan): diagonal nearest-neighbor bonds
   (t\ :sub:`1`, V\ :sub:`1`, J\ :sub:`1` and
   t\ :sub:`1`\ ', V\ :sub:`1`\ ', J\ :sub:`1`\ ').
   Dashed (orange): next-nearest leg bonds
   (t\ :sub:`2`\ ', V\ :sub:`2`\ ', J\ :sub:`2`\ ').

==============  ================================
Property        Value
==============  ================================
Canonical name  ``ladder``
Aliases         ``ladderlattice``
Dimensions      1D
Boost mode      Supported (see constraints below)
==============  ================================

Interactions
------------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Bond class
     - Directions
     - Keywords
   * - Nearest neighbor — rung (W direction)
     - 0
     - ``t0``, ``V0``, ``J0``
   * - Nearest neighbor — diagonal
     - 1
     - ``t1``, ``V1``, ``J1``
   * - Nearest neighbor — leg (L direction)
     - 2
     - ``t2``, ``V2``, ``J2``
   * - 2nd nearest neighbor — diagonal
     - 1
     - ``t1'``, ``V1'``, ``J1'``
   * - 2nd nearest neighbor — leg
     - 2
     - ``t2'``, ``V2'``, ``J2'``
   * - On-site
     - —
     - ``U``, ``mu``, ``h``, ``Gamma``, ``Gamma_y``, ``D``

The ladder lattice does not support 3rd nearest-neighbor interactions.

Boost Mode Constraints
----------------------

- ``S2 = 1`` (S = 1/2 only)
- ``W = 2``
- ``L`` must be even and ≥ 4
