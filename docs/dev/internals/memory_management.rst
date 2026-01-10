Memory Management
=================

This document describes how StdFace allocates and deallocates dynamic memory.

Overview
--------

StdFace uses the ``setmemory.h`` module for all dynamic memory allocation. This
module provides type-safe wrapper functions around standard C ``calloc()`` and
``free()`` calls.

**Source reference**: ``src/setmemory.h``, ``src/setmemory.c``

Design Principles
-----------------

1. **Type safety**: Each data type has dedicated allocation/deallocation functions
2. **Zero initialization**: All allocated memory is zero-initialized via ``calloc()``
3. **Contiguous 2D/3D arrays**: Multi-dimensional arrays use a single contiguous
   memory block with pointer indexing for cache efficiency
4. **Symmetric naming**: Each ``*_allocate()`` function has a matching ``free_*_allocate()``

Naming Convention
-----------------

Function names follow the pattern::

    <type>_<dimension>d_allocate
    free_<type>_<dimension>d_allocate

Where:

- ``<type>``: Data type abbreviation

  - ``ui`` = ``unsigned int``
  - ``li`` = ``long int``
  - ``lui`` = ``long unsigned int``
  - ``i`` = ``int``
  - ``d`` = ``double``
  - ``cd`` = ``complex double``

- ``<dimension>``: Array dimensionality (1, 2, or 3)

Allocation Functions
--------------------

1D Arrays
^^^^^^^^^

.. code-block:: c

   unsigned int     *ui_1d_allocate(const long unsigned int N);
   long int         *li_1d_allocate(const long unsigned int N);
   long unsigned int *lui_1d_allocate(const long unsigned int N);
   int              *i_1d_allocate(const long unsigned int N);
   double           *d_1d_allocate(const long unsigned int N);
   complex double   *cd_1d_allocate(const long unsigned int N);

**Example implementation** (``src/setmemory.c`` lines 112-117):

.. code-block:: c

   int *i_1d_allocate(const long unsigned int N){
       int *A;
       A = (int*)calloc((N), sizeof(int));
       return A;
   }

2D Arrays
^^^^^^^^^

.. code-block:: c

   long int       **li_2d_allocate(const long unsigned int N, const long unsigned int M);
   int            **i_2d_allocate(const long unsigned int N, const long unsigned int M);
   double         **d_2d_allocate(const long unsigned int N, const long unsigned int M);
   complex double **cd_2d_allocate(const long unsigned int N, const long unsigned int M);

**Implementation pattern** (``src/setmemory.c`` lines 133-143):

.. code-block:: c

   int **i_2d_allocate(const long unsigned int N, const long unsigned int M) {
       int **A;
       long unsigned int int_i;
       A = (int **) calloc((N), sizeof(int *));
       A[0] = (int *) calloc((M * N), sizeof(int));
       for (int_i = 0; int_i < N; int_i++) {
           A[int_i] = A[0] + int_i * M;
       }
       return A;
   }

**Memory layout**: The 2D array is stored as a single contiguous block. ``A[0]``
points to the entire data block, and ``A[i]`` pointers are set up for row access.

3D Arrays
^^^^^^^^^

.. code-block:: c

   int            ***i_3d_allocate(const long unsigned int N, const long unsigned int M,
                                   const long unsigned int L);
   complex double ***cd_3d_allocate(const long unsigned int N, const long unsigned int M,
                                    const long unsigned int L);

Deallocation Functions
----------------------

Each allocation function has a corresponding deallocation function:

.. code-block:: c

   void free_ui_1d_allocate(unsigned int *A);
   void free_li_1d_allocate(long int *A);
   void free_li_2d_allocate(long int **A);
   void free_lui_1d_allocate(long unsigned int *A);
   void free_i_1d_allocate(int *A);
   void free_i_2d_allocate(int **A);
   void free_i_3d_allocate(int ***A);
   void free_d_1d_allocate(double *A);
   void free_d_2d_allocate(double **A);
   void free_cd_1d_allocate(double complex *A);
   void free_cd_2d_allocate(double complex **A);
   void free_cd_3d_allocate(double complex ***A);

**Example** (``src/setmemory.c`` lines 122-124):

.. code-block:: c

   void free_i_1d_allocate(int *A){
       free(A);
   }

For 2D and 3D arrays, the deallocation functions free the data block and pointer
arrays in the correct order.

Function Summary
----------------

=======================  ===========  =======================
Allocation Function      Dimension    Deallocation Function
=======================  ===========  =======================
``ui_1d_allocate``       1D           ``free_ui_1d_allocate``
``li_1d_allocate``       1D           ``free_li_1d_allocate``
``li_2d_allocate``       2D           ``free_li_2d_allocate``
``lui_1d_allocate``      1D           ``free_lui_1d_allocate``
``i_1d_allocate``        1D           ``free_i_1d_allocate``
``i_2d_allocate``        2D           ``free_i_2d_allocate``
``i_3d_allocate``        3D           ``free_i_3d_allocate``
``d_1d_allocate``        1D           ``free_d_1d_allocate``
``d_2d_allocate``        2D           ``free_d_2d_allocate``
``cd_1d_allocate``       1D           ``free_cd_1d_allocate``
``cd_2d_allocate``       2D           ``free_cd_2d_allocate``
``cd_3d_allocate``       3D           ``free_cd_3d_allocate``
=======================  ===========  =======================

Usage in StdFace
----------------

The ``StdFace_MallocInteractions()`` function uses these allocation functions
to set up the main data arrays:

**Source reference**: ``src/StdFace_ModelUtil.c``

- ``transindx`` (transfer indices): allocated as 2D ``int`` array
- ``trans`` (transfer coefficients): allocated as 1D ``complex double`` array
- ``intrindx`` (interaction indices): allocated as 2D ``int`` array
- ``intr`` (interaction coefficients): allocated as 1D ``complex double`` array

Ownership and Lifetime
----------------------

- Allocation occurs in ``StdFace_MallocInteractions()`` called from lattice
  constructors
- The ``StdIntList`` structure holds pointers to allocated arrays
- Deallocation occurs at the end of ``StdFace_main()`` after output generation
- Caller (``StdFace_main()``) is responsible for freeing all allocated memory

Error Handling
--------------

The allocation functions do not check for ``NULL`` returns from ``calloc()``.
If memory allocation fails, the behavior is undefined (likely a crash on first
access).

This is an observed pattern; explicit error handling for allocation failure is
not implemented in the current codebase.

Source References
-----------------

- Header: ``src/setmemory.h``
- Implementation: ``src/setmemory.c``
- Usage: ``StdFace_MallocInteractions()`` in ``src/StdFace_ModelUtil.c``
