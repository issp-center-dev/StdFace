Lattice Reference
=================

StdFace supports the following built-in lattice types.

.. _lattice-types:

.. list-table::
   :header-rows: 1
   :widths: 20 12 55 13

   * - Lattice
     - Dimensions
     - Description
     - Boost
   * - :doc:`chain`
     - 1D
     - Linear chain
     - ✓
   * - :doc:`ladder`
     - 1D
     - Two-leg or multi-leg ladder
     - ✓
   * - :doc:`square`
     - 2D
     - Square lattice
     -
   * - :doc:`triangular`
     - 2D
     - Triangular lattice
     -
   * - :doc:`honeycomb`
     - 2D
     - Honeycomb (hexagonal) lattice
     - ✓
   * - :doc:`kagome`
     - 2D
     - Kagome lattice
     - ✓
   * - :doc:`orthorhombic`
     - 3D
     - Orthorhombic / cubic lattice
     -
   * - :doc:`fcortho`
     - 3D
     - Face-centered orthorhombic lattice
     -
   * - :doc:`pyrochlore`
     - 3D
     - Pyrochlore lattice
     -
   * - :doc:`wannier90`
     - Any
     - Import from Wannier90 format files
     -

✓ = HPhi Boost mode supported.

.. toctree::
   :hidden:

   chain
   ladder
   square
   triangular
   honeycomb
   kagome
   orthorhombic
   fcortho
   pyrochlore
   wannier90
