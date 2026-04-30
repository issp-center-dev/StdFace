HPhi Parameters
===============

.. contents::
   :local:
   :depth: 2

Parameters specific to the HPhi solver (``hphi_dry.out``).  For parameters
shared across solvers see :doc:`common`.

.. _method-hphi:

``method``
----------

Specifies the diagonalisation algorithm.

====================  =========================================================
Value                 Algorithm
====================  =========================================================
``"lanczos"``         Lanczos method
``"lanczosenergy"``   Lanczos with eigenvector output
``"tpq"``             Thermal Pure Quantum states
``"fulldiag"``        Full (exact) diagonalisation (aliases: ``"direct"``,
                      ``"alldiag"``)
``"cg"``              Conjugate Gradient
``"timeevolution"``   Real-time evolution (aliases: ``"te"``,
                      ``"time-evolution"``)
``"ctpq"``            Canonical TPQ
====================  =========================================================

``exct``
--------

Number of eigenvalues/eigenvectors to compute (ground state + excited states).
Default: ``1`` (ground state only).

Boost Mode
----------

Boost mode is a high-performance algorithm for S=1/2 spin models that
exploits a 6-spin interaction structure to accelerate the Lanczos
diagonalisation.  Internally, HPhi reorganises the Hilbert space around
*pivot sites* and applies quantum-number projection, achieving better cache
efficiency than the standard algorithm on supported lattices.

Enabling Boost Mode
^^^^^^^^^^^^^^^^^^^

Set ``model`` to one of:

- ``"spingcboost"`` — grand-canonical spin with Boost
- ``"spingccma"`` — grand-canonical spin with CMA variant

.. code-block:: text

   model   = "spingcboost"
   lattice = chain
   L       = 16

Supported Lattices and Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+---------------+--------------------------------------------------+
| Lattice       | Constraints                                      |
+===============+==================================================+
| ``chain``     | ``S2 = 1`` (S=1/2); ``L`` must be a multiple of 8|
+---------------+--------------------------------------------------+
| ``ladder``    | ``S2 = 1``; ``W = 2``; ``L`` even and ≥ 4        |
+---------------+--------------------------------------------------+
| ``honeycomb`` | ``S2 = 1``                                       |
+---------------+--------------------------------------------------+
| ``kagome``    | ``S2 = 1``                                       |
+---------------+--------------------------------------------------+

Other lattices do not support Boost mode and will produce an error.

Additional Output File
^^^^^^^^^^^^^^^^^^^^^^

When Boost mode is active, ``boost.def`` is generated in addition to the
standard HPhi files.  It contains the magnetic field values, exchange
interaction matrices, and the pivot-site topology (pivot count, 6-spin
pair lists, spin-shift information) required by HPhi's Boost kernel.
