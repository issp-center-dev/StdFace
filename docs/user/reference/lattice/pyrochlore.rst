Pyrochlore Lattice
==================

Three-dimensional pyrochlore lattice.  The unit cell contains four sites
arranged in corner-sharing tetrahedra.

==============  =====================================
Property        Value
==============  =====================================
Canonical name  ``pyrochlore``
Aliases         —
Dimensions      3D
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
   * - Nearest neighbor
     - 0, 1, 2
     - ``t0``, ``V0``, ``J0``; ``t1``, ``V1``, ``J1``; ``t2``, ``V2``, ``J2``
   * - 2nd nearest neighbor
     - 0, 1, 2
     - ``t0'``, ``V0'``, ``J0'``; ``t1'``, ``V1'``, ``J1'``; ``t2'``, ``V2'``, ``J2'``
   * - On-site
     - —
     - ``U``, ``mu``, ``h``, ``Gamma``, ``Gamma_y``, ``D``

The pyrochlore lattice does not support 3rd nearest-neighbor interactions.
