Testing
=======

StdFace has two complementary test suites: integration tests that exercise
the compiled C binaries against reference output, and unit tests that cover
the Python reimplementation.

.. contents::
   :local:
   :depth: 2

Test Types
----------

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Type
     - Location
     - Framework
   * - Integration tests
     - ``test/hphi/``, ``test/mvmc/``, ``test/uhf/``, ``test/hwave/``
     - CMake / CTest + shell scripts
   * - Unit tests
     - ``test/unit/``
     - pytest (Python implementation only)

Integration Tests
-----------------

Running
^^^^^^^

Integration tests are registered with CTest during the CMake configuration
step.  Enable the solver(s) you want to test:

.. code-block:: bash

   mkdir build && cd build
   cmake .. -DHPHI=ON -DMVMC=ON -DUHF=ON -DHWAVE=ON
   make
   ctest

Useful CTest options:

.. code-block:: bash

   ctest -R "^hphi_" -V      # HPhi tests only, verbose output
   ctest -R "^mvmc_" -V      # mVMC tests only
   ctest -j4                  # Run 4 tests in parallel

Test Case Structure
^^^^^^^^^^^^^^^^^^^

Each test case is a directory containing an input file and a ``ref/``
subdirectory of expected output files:

.. code-block:: text

   test/SOLVER/TEST_NAME/
   ├── stan.in        # StdFace input file (StdFace.def for mVMC)
   └── ref/           # Reference output files
       ├── namelist.def
       ├── modpara.def
       ├── trans.def
       ├── lattice.gp
       └── ...

The test harness (``test/common/check_case.sh``) copies the test directory
to a temporary working directory, runs the dry-run generator, and compares
every file in ``ref/`` against the generated output.

Output Comparison and Tolerance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Exact text files are compared with ``diff``.  Files containing floating-point
numbers are compared with ``test/common/almost_diff.py``, which accepts a
relative error up to **1.0 × 10⁻⁸**.  This handles platform-specific
variations in the least-significant digits of floating-point output.
``NaN`` vs ``-NaN`` differences are also tolerated.

TEST_MODE
^^^^^^^^^

The ``TEST_MODE`` variable in ``test/CMakeLists.txt`` controls which
generator is under test:

==============  =======================================================
Value           Generator
==============  =======================================================
``"base"``      Compiled C binaries (``hphi_dry.out``, etc.) — default
``"python"``    Python CLI (``python3 -mstdface --solver SOLVER``)
==============  =======================================================

Change the value before running ``cmake`` to switch between the two
implementations.

Test Categories
^^^^^^^^^^^^^^^

Tests are grouped by name prefix:

=====================  =====================================================
Prefix / Pattern       Category
=====================  =====================================================
``matrix_*``           Auto-generated cases covering all geometry/model
                       combinations (e.g. ``matrix_chain_hubbard``)
``fail_*``             Negative tests — verify that invalid input causes
                       a non-zero exit status (e.g. ``fail_unsupported_keyword``,
                       ``fail_kagome_boost_jprime``)
``boost_*``            HPhi Boost-mode cases
``ctpq_``, ``tpq_``    Finite-temperature TPQ cases
``fulldiag_``          Full diagonalisation cases
``te_*``               Time-evolution cases
``spectrum_*``         Spectral function cases
``boundary_*``         Boundary condition cases
=====================  =====================================================

Regenerating Reference Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the expected output of StdFace changes (e.g. after a bug fix or new
feature), regenerate the reference files with the solver-specific script:

.. code-block:: bash

   test/hphi/make_reference.sh
   test/mvmc/make_reference.sh
   test/uhf/make_reference.sh
   test/hwave/make_reference.sh

These scripts run each test case and overwrite the contents of the
corresponding ``ref/`` directory.  Commit the updated reference files
together with the source change.

Unit Tests
----------

The unit tests cover the Python reimplementation of StdFace
(``python/stdface/``).  They do not require a compiled binary.

Running
^^^^^^^

.. code-block:: bash

   pytest test/unit/

.. code-block:: bash

   pytest test/unit/ -v          # verbose
   pytest test/unit/ -k lattice  # run tests matching "lattice"

Structure
^^^^^^^^^

.. code-block:: text

   test/unit/
   ├── conftest.py               # adds python/ to sys.path
   ├── test_stdface_vals.py      # core data structures
   ├── test_keyword_parser.py    # input file parsing
   ├── test_interaction_builder.py
   ├── test_lattices.py
   ├── test_hphi_writer.py
   ├── test_mvmc_writer.py
   ├── test_common_writer.py
   ├── test_wannier90.py
   └── ...                       # 28 files total

Test counts: approximately 366 integration test cases and 28 unit test
files across all solvers.
