Error Handling
==============

This document describes how StdFace handles errors and abnormal conditions.

Overview
--------

StdFace uses a simple error handling strategy:

1. Detect invalid or missing parameters during input processing
2. Print a descriptive error message to ``stdout``
3. Terminate the program via ``StdFace_exit()``

There are no exceptions or error recovery mechanisms. All errors are fatal and
result in program termination.

Exit Function
-------------

The central error exit mechanism is ``StdFace_exit()``:

**Source reference**: ``src/StdFace_ModelUtil.c`` lines 34-47

.. code-block:: c

   void StdFace_exit(int errorcode)
   {
     int ierr = 0;
     fflush(stdout);
     fflush(stderr);
   #ifdef MPI
     fprintf(stdout, "\n\n #######  You DO NOT have to WORRY about "
                     "the following MPI-ERROR MESSAGE.  #######\n\n");
     ierr = MPI_Abort(MPI_COMM_WORLD, errorcode);
     ierr = MPI_Finalize();
   #endif
     if (ierr != 0) fprintf(stderr, "\n  MPI_Finalize() = %d\n\n", ierr);
     exit(errorcode);
   }

**Behavior**:

- Flushes ``stdout`` and ``stderr`` to ensure all messages are written
- If compiled with MPI support (``#ifdef MPI``):

  - Prints a reassurance message about MPI error output
  - Calls ``MPI_Abort()`` to terminate all MPI processes
  - Calls ``MPI_Finalize()`` for cleanup

- Terminates the program with the specified ``errorcode``

Error Categories
----------------

Required Parameter Missing
^^^^^^^^^^^^^^^^^^^^^^^^^^

When a required input parameter is not specified, ``StdFace_RequiredVal_i()``
reports the error and exits.

**Source reference**: ``src/StdFace_ModelUtil.c`` lines 588-600

.. code-block:: c

   void StdFace_RequiredVal_i(char* valname, int val)
   {
     int NaN_i = 2147483647;

     if (val == NaN_i){
       fprintf(stdout, "ERROR ! %s is NOT specified !\n", valname);
       StdFace_exit(-1);
     }
     else fprintf(stdout, "  %15s = %-3d\n", valname, val);
   }

**Detection mechanism**: Integer parameters are initialized to ``2147483647``
(``INT_MAX``) as a sentinel value. If the value remains unchanged after input
parsing, the parameter was not specified.

Unused Parameter Specified
^^^^^^^^^^^^^^^^^^^^^^^^^^

When a user specifies a parameter that is not applicable to the current model
configuration, the ``StdFace_NotUsed_*`` family of functions reports this as an
error.

**Source reference**: ``src/StdFace_ModelUtil.c`` lines 503-514

.. code-block:: c

   void StdFace_NotUsed_d(char* valname, double val)
   {
     if (isnan(val) == 0) {
       fprintf(stdout, "\n Check !  %s is SPECIFIED but will NOT be USED. \n",
               valname);
       fprintf(stdout, "            Please COMMENT-OUT this line \n");
       fprintf(stdout, "            or check this input is REALLY APPROPRIATE "
                       "for your purpose ! \n\n");
       StdFace_exit(-1);
     }
   }

**Variants**:

- ``StdFace_NotUsed_d()`` - for ``double`` parameters (uses ``isnan()`` check)
- ``StdFace_NotUsed_i()`` - for ``int`` parameters (uses sentinel check)
- ``StdFace_NotUsed_c()`` - for ``complex double`` parameters
- ``StdFace_NotUsed_J()`` - for exchange coupling parameters (J matrices)

**Detection mechanism**: Double parameters are initialized to ``NaN``. If the
value is not ``NaN`` after input parsing, the parameter was specified.

Exit Codes
----------

The following exit codes are used:

=============  =====================================================
Code           Meaning
=============  =====================================================
``-1``         Required parameter missing or unused parameter error
``errorcode``  Passed to ``StdFace_exit()`` from caller
=============  =====================================================

Error Output Format
-------------------

Error messages are written to ``stdout`` (not ``stderr``). Observed patterns:

- Required parameter missing: ``ERROR ! <name> is NOT specified !``
- Unused parameter: Multi-line message suggesting to comment out the parameter

Design Rationale
----------------

The error handling approach prioritizes:

1. **Simplicity**: No exception handling overhead or complex error propagation
2. **User guidance**: Error messages tell users how to fix the problem
3. **Fail-fast**: Invalid configurations are caught early during input parsing
4. **MPI awareness**: Proper cleanup of MPI environment on error

Source References
-----------------

- Exit function: ``src/StdFace_ModelUtil.c`` lines 34-47
- Required value check: ``src/StdFace_ModelUtil.c`` lines 588-600
- Not-used checks: ``src/StdFace_ModelUtil.c`` lines 503-569
- Declaration: ``src/StdFace_ModelUtil.h``
