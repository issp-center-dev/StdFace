Execution Flow
==============

High-Level Flow
---------------

The execution flow of StdFace follows a linear pipeline:

::

   CLI Invocation
        │
        ▼
   ┌─────────────────────┐
   │ main() in dry.c     │ ─── Argument parsing
   └─────────┬───────────┘
             │
             ▼
   ┌─────────────────────┐
   │ StdFace_main()      │ ─── Core processing
   │ in StdFace_main.c   │
   └─────────┬───────────┘
             │
             ├──► Input file parsing
             │
             ├──► Lattice construction
             │
             └──► Output file generation
                        │
                        ▼
                  Configuration files
                  for target solver

**Source reference**: ``src/dry.c``, ``src/StdFace_main.c``

CLI Entry Point
---------------

The ``main()`` function in ``src/dry.c`` handles:

.. code-block:: c

   int main(int argc, char *argv[])
   {
     int status = 0;
     if (argc < 2) {
       printVersion();           // Display version
       usage(argv[0]);           // Show usage
       status = 1;
     } else if (strcmp(argv[1], "-v") == 0) {
       printVersion();           // Version flag
     } else {
       StdFace_main(argv[1]);    // Process input file
     }
     return status;
   }

**Behavior**:

1. No arguments: Print version and usage, exit with status 1
2. ``-v`` flag: Print version only
3. Input file provided: Dispatch to ``StdFace_main()``

**Source reference**: ``src/dry.c:34-47``

StdFace_main() Processing
-------------------------

The ``StdFace_main()`` function (lines 2454-3066 of ``src/StdFace_main.c``)
performs the following steps:

1. Initialization
^^^^^^^^^^^^^^^^^

.. code-block:: c

   struct StdIntList *StdI;
   StdI = (struct StdIntList *)malloc(sizeof(struct StdIntList));
   StdFace_ResetVals(StdI);

- Allocates the central ``StdIntList`` structure
- Initializes all fields to default/sentinel values via ``StdFace_ResetVals()``

**Source reference**: ``src/StdFace_main.c:2459-2471``

2. Input File Parsing
^^^^^^^^^^^^^^^^^^^^^

The input file is parsed line by line:

.. code-block:: c

   while (fgets(ctmpline, 256, fp) != NULL) {
     TrimSpaceQuote(ctmpline);
     // Skip comments and empty lines
     if (strncmp(ctmpline, "//", 2) == 0) continue;
     if (ctmpline[0] == '\0') continue;

     keyword = strtok(ctmpline, "=");
     value = strtok(NULL, "=");
     Text2Lower(keyword);

     // Match keyword and store value
     if (strcmp(keyword, "lattice") == 0) StoreWithCheckDup_sl(...);
     else if (strcmp(keyword, "model") == 0) StoreWithCheckDup_sl(...);
     // ... many more keyword handlers
   }

**Key parsing functions**:

+--------------------------+------------------------------------------+
| Function                 | Purpose                                  |
+==========================+==========================================+
| ``TrimSpaceQuote()``     | Remove whitespace and quotes from line   |
+--------------------------+------------------------------------------+
| ``Text2Lower()``         | Convert keyword to lowercase             |
+--------------------------+------------------------------------------+
| ``StoreWithCheckDup_i``  | Store integer with duplicate check       |
+--------------------------+------------------------------------------+
| ``StoreWithCheckDup_d``  | Store double with duplicate check        |
+--------------------------+------------------------------------------+
| ``StoreWithCheckDup_c``  | Store complex with duplicate check       |
+--------------------------+------------------------------------------+
| ``StoreWithCheckDup_s``  | Store string with duplicate check        |
+--------------------------+------------------------------------------+
| ``StoreWithCheckDup_sl`` | Store string (lowercase) with check      |
+--------------------------+------------------------------------------+

**Source reference**: ``src/StdFace_main.c:2482-2776`` (parsing loop)

3. Model Validation
^^^^^^^^^^^^^^^^^^^

After parsing, the model type is validated and normalized:

.. code-block:: c

   if (strcmp(StdI->model, "fermionhubbard") == 0
     || strcmp(StdI->model, "hubbard") == 0)
     strcpy(StdI->model, "hubbard\0");
   else if (strcmp(StdI->model, "spin") == 0)
     strcpy(StdI->model, "spin\0");
   else if (strcmp(StdI->model, "kondolattice") == 0
     || strcmp(StdI->model, "kondo") == 0)
     strcpy(StdI->model, "kondo\0");
   // Grand canonical variants set StdI->lGC = 1

**Supported models**:

- ``hubbard`` / ``fermionhubbard`` (canonical or grand canonical)
- ``spin`` (canonical or grand canonical)
- ``kondo`` / ``kondolattice`` (canonical or grand canonical)

**Source reference**: ``src/StdFace_main.c:2795-2830``

4. Lattice Construction
^^^^^^^^^^^^^^^^^^^^^^^

Based on the ``lattice`` keyword, the appropriate constructor is called:

.. code-block:: c

   if (strcmp(StdI->lattice, "chain") == 0)
     StdFace_Chain(StdI);
   else if (strcmp(StdI->lattice, "square") == 0
     || strcmp(StdI->lattice, "tetragonal") == 0)
     StdFace_Tetragonal(StdI);
   else if (strcmp(StdI->lattice, "honeycomb") == 0)
     StdFace_Honeycomb(StdI);
   // ... other lattice types

Each lattice constructor:

1. Sets up site geometry
2. Defines neighbor relationships
3. Computes hopping and interaction terms
4. Populates ``StdI->trans[]`` and ``StdI->intr[]`` arrays

**Source reference**: ``src/StdFace_main.c:2915-2945``

5. Mode-Specific Output Generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Output generation differs based on the compile-time mode:

**HPhi mode** (``_HPhi`` defined):

.. code-block:: c

   #if defined(_HPhi)
   PrintLocSpin(StdI);
   PrintTrans(StdI);
   PrintInteractions(StdI);
   CheckModPara(StdI);
   PrintModPara(StdI);
   PrintExcitation(StdI);
   PrintCalcMod(StdI);
   CheckOutputMode(StdI);
   Print1Green(StdI);
   Print2Green(StdI);
   PrintNamelist(StdI);
   #endif

**mVMC mode** (``_mVMC`` defined):

- Similar structure but includes variational parameters
- Generates Jastrow, Gutzwiller, orbital files

**UHF mode** (``_UHF`` defined):

- Simpler output; no Green's function files

**HWAVE mode** (``_HWAVE`` defined):

- Two sub-modes: ``uhfr`` (real-space) or ``uhfk``/``rpa`` (k-space)
- K-space modes call ``ExportGeometry()`` and ``ExportInteraction()``

**Source reference**: ``src/StdFace_main.c:2960-3040``

6. Cleanup
^^^^^^^^^^

Memory allocated during processing is freed:

.. code-block:: c

   free(StdI->locspinflag);
   for (ktrans = 0; ktrans < StdI->ntrans; ktrans++) {
     free(StdI->transindx[ktrans]);
   }
   free(StdI->transindx);
   free(StdI->trans);
   // ... similar for other arrays
   free(StdI);

**Source reference**: ``src/StdFace_main.c:3050-3065``

Mode-Specific Behavior
----------------------

The same source code behaves differently based on compile-time defines:

+-------------+----------------------------------------------------------+
| Mode        | Key Differences                                          |
+=============+==========================================================+
| ``_HPhi``   | Full output including spectrum calculation parameters,   |
|             | TPQ settings, Green's functions                          |
+-------------+----------------------------------------------------------+
| ``_mVMC``   | Includes variational parameters (Jastrow, Gutzwiller),   |
|             | sublattice settings for symmetry                         |
+-------------+----------------------------------------------------------+
| ``_UHF``    | Mean-field settings, iteration parameters, mixing        |
+-------------+----------------------------------------------------------+
| ``_HWAVE``  | Two sub-modes: real-space (like UHF) or k-space          |
|             | (Wannier90-style export)                                 |
+-------------+----------------------------------------------------------+

The mode affects:

1. Which keywords are valid in the input file
2. What output files are generated
3. Which calculation parameters are required

**Source reference**: ``src/StdFace_main.c:2664-2774`` (mode-specific keywords)

Error Handling
--------------

Errors are handled by the ``StdFace_exit()`` function:

- Prints error message to stdout
- Calls ``exit()`` with the provided status code

Common error conditions:

- Missing required parameters
- Invalid keyword values
- Unsupported model/lattice combinations
- File open failures

Error messages are printed to stdout (not stderr). This is the observed pattern
in the codebase.

**Source reference**: ``src/StdFace_ModelUtil.h`` (``StdFace_exit``),
``src/StdFace_main.c`` (usage examples)
