Validation and Limitations
==========================

.. contents::
   :local:
   :depth: 2

Validation and Errors
---------------------

StdFace performs input validation during processing and reports errors to
``stdout``.

Common Error Messages
^^^^^^^^^^^^^^^^^^^^^

**ERROR ! <parameter> is NOT specified !**

A parameter expected by StdFace for the current model/lattice configuration
was not found in the input file. Add the indicated parameter to your input.

*Cause*: Integer parameters are initialized to a sentinel value (2147483647).
If unchanged after parsing, the parameter is considered missing.

**ERROR ! Keyword <parameter> is duplicated !**

The same keyword appears more than once in the input file.  StdFace does not
accept duplicate keys — it exits immediately with status ``-1``.

*Resolution*: Remove the duplicate entry, keeping only one occurrence.

**Check ! <parameter> is SPECIFIED but will NOT be USED.**

A parameter was specified in the input file that does not apply to the current
model/lattice combination.

*Resolution*: Remove or comment out the unused parameter from your input file.

*Cause*: Double parameters are initialized to NaN. If a non-NaN value is
detected for an inapplicable parameter, this error is raised.

What StdFace Validates
^^^^^^^^^^^^^^^^^^^^^^

- Presence of expected parameters (via sentinel value detection)
- Absence of inapplicable parameters (via NaN detection for doubles)
- Uniqueness of each keyword (duplicate keys are a fatal error)

What StdFace Does NOT Validate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Physical consistency of parameter values
- Numerical stability or correctness of the model
- Solver-specific parameter ranges or constraints
- Completeness of Wannier90 input files

Limitations and Non-Goals
--------------------------

StdFace is an input file generator, not a solver. The following are outside
its scope.

Solver-Specific Semantics
^^^^^^^^^^^^^^^^^^^^^^^^^

StdFace generates input files but does not interpret solver-specific options
beyond what is needed for file generation. Parameters like ``method`` and
``exct`` are passed through to output files; their validity is not checked
by StdFace.

Numerical Correctness
^^^^^^^^^^^^^^^^^^^^^

StdFace does not validate that the specified model parameters produce
physically meaningful or numerically stable results. Users are responsible
for ensuring appropriate parameter values.

Advanced Solver Options
^^^^^^^^^^^^^^^^^^^^^^^

Some solver options may not be exposed through StdFace's input format. For
advanced configuration, users may need to edit the generated files directly
or consult solver-specific documentation.

Output File Formats
^^^^^^^^^^^^^^^^^^^

The format and content of output files follow the input specifications of
each target solver.  For details of each file's internal format, refer to
the documentation of HPhi, mVMC, UHF, or H-wave respectively.
