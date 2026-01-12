setmemory.h
===========

Purpose
-------

Memory allocation and deallocation utilities for various array types. Provides
type-safe allocation functions for 1D, 2D, and 3D arrays of integers, doubles,
and complex doubles.

**Source reference**: ``src/setmemory.h``

Public Types
------------

None.

Public Macros/Constants
-----------------------

.. c:macro:: MVMC_SETMEMORY_H

   Include guard macro.

Public Functions
----------------

Unsigned Integer Arrays
^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: unsigned int *ui_1d_allocate(const long unsigned int N)

   Allocates a 1D array of unsigned integers.

   :param N: Size of the array
   :returns: Pointer to allocated array

.. c:function:: void free_ui_1d_allocate(unsigned int *A)

   Frees a 1D unsigned integer array.

   :param A: Pointer to array to free

Long Integer Arrays
^^^^^^^^^^^^^^^^^^^

.. c:function:: long int *li_1d_allocate(const long unsigned int N)

   Allocates a 1D array of long integers.

   :param N: Size of the array
   :returns: Pointer to allocated array

.. c:function:: void free_li_1d_allocate(long int *A)

   Frees a 1D long integer array.

   :param A: Pointer to array to free

.. c:function:: long int **li_2d_allocate(const long unsigned int N, const long unsigned int M)

   Allocates a 2D array of long integers.

   :param N: First dimension size
   :param M: Second dimension size
   :returns: Pointer to allocated 2D array

.. c:function:: void free_li_2d_allocate(long int **A)

   Frees a 2D long integer array.

   :param A: Pointer to 2D array to free

Long Unsigned Integer Arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: long unsigned int *lui_1d_allocate(const long unsigned int N)

   Allocates a 1D array of long unsigned integers.

   :param N: Size of the array
   :returns: Pointer to allocated array

.. c:function:: void free_lui_1d_allocate(long unsigned int *A)

   Frees a 1D long unsigned integer array.

   :param A: Pointer to array to free

Integer Arrays
^^^^^^^^^^^^^^

.. c:function:: int *i_1d_allocate(const long unsigned int N)

   Allocates a 1D array of integers.

   :param N: Size of the array
   :returns: Pointer to allocated array

.. c:function:: void free_i_1d_allocate(int *A)

   Frees a 1D integer array.

   :param A: Pointer to array to free

.. c:function:: int **i_2d_allocate(const long unsigned int N, const long unsigned int M)

   Allocates a 2D array of integers.

   :param N: First dimension size
   :param M: Second dimension size
   :returns: Pointer to allocated 2D array

.. c:function:: void free_i_2d_allocate(int **A)

   Frees a 2D integer array.

   :param A: Pointer to 2D array to free

.. c:function:: int ***i_3d_allocate(const long unsigned int N, const long unsigned int M, const long unsigned int L)

   Allocates a 3D array of integers.

   :param N: First dimension size
   :param M: Second dimension size
   :param L: Third dimension size
   :returns: Pointer to allocated 3D array

.. c:function:: void free_i_3d_allocate(int ***A)

   Frees a 3D integer array.

   :param A: Pointer to 3D array to free

Double Arrays
^^^^^^^^^^^^^

.. c:function:: double *d_1d_allocate(const long unsigned int N)

   Allocates a 1D array of doubles.

   :param N: Size of the array
   :returns: Pointer to allocated array

.. c:function:: void free_d_1d_allocate(double *A)

   Frees a 1D double array.

   :param A: Pointer to array to free

.. c:function:: double **d_2d_allocate(const long unsigned int N, const long unsigned int M)

   Allocates a 2D array of doubles.

   :param N: First dimension size
   :param M: Second dimension size
   :returns: Pointer to allocated 2D array

.. c:function:: void free_d_2d_allocate(double **A)

   Frees a 2D double array.

   :param A: Pointer to 2D array to free

Complex Double Arrays
^^^^^^^^^^^^^^^^^^^^^

.. c:function:: complex double *cd_1d_allocate(const long unsigned int N)

   Allocates a 1D array of complex doubles.

   :param N: Size of the array
   :returns: Pointer to allocated array

.. c:function:: void free_cd_1d_allocate(double complex *A)

   Frees a 1D complex double array.

   :param A: Pointer to array to free

.. c:function:: complex double **cd_2d_allocate(const long unsigned int N, const long unsigned int M)

   Allocates a 2D array of complex doubles.

   :param N: First dimension size
   :param M: Second dimension size
   :returns: Pointer to allocated 2D array

.. c:function:: void free_cd_2d_allocate(double complex **A)

   Frees a 2D complex double array.

   :param A: Pointer to 2D array to free

.. c:function:: double complex ***cd_3d_allocate(const long unsigned int N, const long unsigned int M, const long unsigned int L)

   Allocates a 3D array of complex doubles.

   :param N: First dimension size
   :param M: Second dimension size
   :param L: Third dimension size
   :returns: Pointer to allocated 3D array

.. c:function:: void free_cd_3d_allocate(double complex ***A)

   Frees a 3D complex double array.

   :param A: Pointer to 3D array to free

Ownership and Lifetime Rules
----------------------------

- Allocation functions return newly allocated memory; caller is responsible for freeing.
- Each ``*_allocate()`` function has a corresponding ``free_*_allocate()`` function.
- Memory is allocated using standard ``malloc()``/``calloc()``.
- Passing NULL to free functions is unspecified in the current code.

Error Handling
--------------

Unspecified in the current code. Allocation failure behavior (e.g., NULL check)
is not documented in the header.

Thread / MPI Safety
-------------------

Unspecified in the current code. Uses standard C memory allocation which is
typically thread-safe in modern implementations.

Source Reference
----------------

- Header: ``src/setmemory.h``
- Implementation: ``src/setmemory.c``
- Author: Kazuyoshi Yoshimi (University of Tokyo)
