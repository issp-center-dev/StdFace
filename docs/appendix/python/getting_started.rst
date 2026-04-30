Getting Started (Python)
========================

.. note::

   This section describes the **Python reimplementation** of StdFace
   (``python/stdface/``).  The main manual covers the C implementation.

Requirements
------------

- Python 3.10 or later
- NumPy

Installation
------------

.. code-block:: bash

   cd python
   pip install -e .           # install with runtime dependencies (numpy)
   pip install -e ".[dev]"    # also install dev dependencies (pytest, pytest-cov)

The editable install (``-e``) is recommended for development; source changes
take effect immediately without reinstalling.

Basic Usage
-----------

Command Line
^^^^^^^^^^^^

After installation, the ``stdface`` command is available:

.. code-block:: bash

   stdface stan.in                  # default solver (HPhi)
   stdface stan.in --solver HPhi
   stdface stan.in --solver mVMC
   stdface stan.in --solver UHF
   stdface stan.in --solver HWAVE
   stdface -v                       # print version

The input file format is identical to the C version (see :doc:`/user/reference/input/common`).

Python API
^^^^^^^^^^

.. code-block:: python

   import sys
   sys.path.insert(0, "python")

   from stdface.core.stdface_main import stdface_main

   stdface_main("stan.in", solver="HPhi")

Running Tests
-------------

.. code-block:: bash

   # Unit tests
   python3 -m pytest test/unit/ -v

   # Unit tests with coverage report
   python3 -m pytest test/unit/ --cov=python --cov-report=html

   # Integration tests (all solvers, all lattices)
   bash test/run_all_integration.sh
