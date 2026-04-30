Chain Lattice
=============

One-dimensional linear chain.

.. figure:: /_static/lattice/chain.png
   :width: 300px
   :align: center

   Chain lattice.  Solid (magenta): nearest-neighbor bonds
   (t\ :sub:`0`, V\ :sub:`0`, J\ :sub:`0`).
   Dashed (orange): next-nearest-neighbor bonds
   (t\ :sub:`0`\ ', V\ :sub:`0`\ ', J\ :sub:`0`\ ').

==============  ================================
Property        Value
==============  ================================
Canonical name  ``chain``
Aliases         ``chainlattice``
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
   * - Nearest neighbor
     - 0
     - ``t0``, ``V0``, ``J0``
   * - 2nd nearest neighbor
     - 0
     - ``t0'``, ``V0'``, ``J0'``
   * - 3rd nearest neighbor
     - 0
     - ``t0''``, ``V0''``, ``J0''``
   * - On-site
     - —
     - ``U``, ``mu``, ``h``, ``Gamma``, ``Gamma_y``, ``D``

Because there is only one bond direction, the isotropic shorthands
``t``, ``V``, ``J`` are equivalent to ``t0``, ``V0``, ``J0``.
Similarly ``t'`` / ``V'`` / ``J'`` and ``t''`` / ``V''`` / ``J''`` are accepted.

Boost Mode Constraints
----------------------

- ``S2 = 1`` (S = 1/2 only)
- ``L`` must be a multiple of 8
