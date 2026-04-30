Honeycomb Lattice
=================

Two-dimensional honeycomb (hexagonal) lattice.  The unit cell contains
two sites (sublattice A and B).

.. figure:: /_static/lattice/honeycomb.png
   :width: 480px
   :align: center

   Anisotropic honeycomb lattice.  Three nearest-neighbor bond directions
   (t\ :sub:`0`, V\ :sub:`0`, J\ :sub:`0`),
   (t\ :sub:`1`, V\ :sub:`1`, J\ :sub:`1`),
   (t\ :sub:`2`, V\ :sub:`2`, J\ :sub:`2`).
   Primed and double-primed labels indicate 2nd and 3rd nearest-neighbor
   bonds respectively.

==============  ================================
Property        Value
==============  ================================
Canonical name  ``honeycomb``
Aliases         ``honeycomblattice``
Dimensions      2D
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
   * - Nearest neighbor
     - 0, 1, 2
     - ``t0``, ``V0``, ``J0``; ``t1``, ``V1``, ``J1``; ``t2``, ``V2``, ``J2``
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

Boost Mode Constraints
----------------------

- ``S2 = 1`` (S = 1/2 only)
