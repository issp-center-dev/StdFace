Orthorhombic Lattice
====================

Three-dimensional orthorhombic (simple cubic) lattice.

.. note::

   In StdFace, ``orthorhombic`` is the general axis-aligned 3D lattice.  The
   alias ``cubic`` is accepted for convenience but does **not** enforce
   W = L = Height — all three dimensions may be set independently.  This
   differs from the standard crystallographic usage where "orthorhombic" and
   "cubic" denote distinct crystal systems.

==============  =====================================
Property        Value
==============  =====================================
Canonical name  ``orthorhombic``
Aliases         ``simplecubic``, ``cubic``
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
   * - Nearest neighbor — along W
     - 0
     - ``t0``, ``V0``, ``J0``
   * - Nearest neighbor — along L
     - 1
     - ``t1``, ``V1``, ``J1``
   * - Nearest neighbor — along H
     - 2
     - ``t2``, ``V2``, ``J2``
   * - 2nd nearest neighbor
     - 0, 1, 2
     - ``t0'``, ``V0'``, ``J0'``; ``t1'``, ``V1'``, ``J1'``; ``t2'``, ``V2'``, ``J2'``
   * - 3rd nearest neighbor (isotropic)
     - —
     - ``t''``, ``V''``, ``J''``
   * - On-site
     - —
     - ``U``, ``mu``, ``h``, ``Gamma``, ``Gamma_y``, ``D``

The 3rd nearest-neighbor interaction accepts only the isotropic scalar form
(``t''``, ``V''``, ``J''``); directional variants ``t0''`` etc. are not
supported for this lattice.
