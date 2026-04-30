Kagome Lattice
==============

Two-dimensional kagome lattice.  The unit cell contains three sites.

.. figure:: /_static/lattice/kagome.png
   :width: 460px
   :align: center

   Kagome lattice.  Three nearest-neighbor bond directions
   (t\ :sub:`0`, V\ :sub:`0`, J\ :sub:`0`),
   (t\ :sub:`1`, V\ :sub:`1`, J\ :sub:`1`),
   (t\ :sub:`2`, V\ :sub:`2`, J\ :sub:`2`).
   Primed labels (dashed) indicate 2nd nearest-neighbor bonds.

==============  ================================
Property        Value
==============  ================================
Canonical name  ``kagome``
Aliases         ``kagomelattice``
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
   * - On-site
     - —
     - ``U``, ``mu``, ``h``, ``Gamma``, ``Gamma_y``, ``D``

The kagome lattice does not support 3rd nearest-neighbor interactions.
The isotropic shorthand ``t`` / ``V`` / ``J`` sets all three nearest-neighbor
directions simultaneously.

Boost Mode Constraints
----------------------

- ``S2 = 1`` (S = 1/2 only)
