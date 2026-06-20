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

The simplest entry point mirrors the command line and writes files to the
current directory:

.. code-block:: python

   from stdface.core.stdface_main import stdface_main

   stdface_main("stan.in", solver="HPhi")

For library use, ``generate()`` returns the output as an object and can write
to a chosen directory (or not write at all):

.. code-block:: python

   from stdface.core.stdface_main import generate

   output = generate("stan.in", solver="HPhi", output_dir="run1")
   data = output.to_dict()        # output_dir=None to skip writing files

Running Tests
-------------

.. code-block:: bash

   # Unit tests
   python3 -m pytest test/unit/ -v

   # Unit tests with coverage report
   python3 -m pytest test/unit/ --cov=python --cov-report=html

   # Integration tests (all solvers, all lattices)
   bash test/run_all_integration.sh
