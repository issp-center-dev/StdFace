/**
 * @file setmemory.h
 * @brief Memory allocation utilities for arrays of various types
 * @copyright 
 * HPhi-mVMC-StdFace - Common input generator
 * Copyright (C) 2015 The University of Tokyo
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MVMC_SETMEMORY_H
#define MVMC_SETMEMORY_H

#include <stdlib.h>
#include <string.h>
#include <complex.h>

/**
 * @brief Allocate 1D array of unsigned integers
 * @param N Size of array to allocate
 * @return Pointer to allocated array, initialized to zero
 */
unsigned int *ui_1d_allocate(const long unsigned int N);

/**
 * @brief Free 1D array of unsigned integers
 * @param A Pointer to array to free
 */
void free_ui_1d_allocate(unsigned int *A);

/**
 * @brief Allocate 1D array of long integers
 * @param N Size of array to allocate
 * @return Pointer to allocated array, initialized to zero
 */
long int *li_1d_allocate(const long unsigned int N);

/**
 * @brief Free 1D array of long integers
 * @param A Pointer to array to free
 */
void free_li_1d_allocate(long int *A);

/**
 * @brief Allocate 2D array of long integers
 * @param N Number of rows
 * @param M Number of columns
 * @return Pointer to allocated array, initialized to zero
 */
long int **li_2d_allocate(const long unsigned int N, const long unsigned int M);

/**
 * @brief Free 2D array of long integers
 * @param A Pointer to array to free
 */
void free_li_2d_allocate(long int **A);

/**
 * @brief Allocate 1D array of long unsigned integers
 * @param N Size of array to allocate
 * @return Pointer to allocated array, initialized to zero
 */
long unsigned int *lui_1d_allocate(const long unsigned int N);

/**
 * @brief Free 1D array of long unsigned integers
 * @param A Pointer to array to free
 */
void free_lui_1d_allocate(long unsigned int *A);

/**
 * @brief Allocate 1D array of integers
 * @param N Size of array to allocate
 * @return Pointer to allocated array, initialized to zero
 */
int *i_1d_allocate(const long unsigned int N);

/**
 * @brief Free 1D array of integers
 * @param A Pointer to array to free
 */
void free_i_1d_allocate(int *A);

/**
 * @brief Allocate 2D array of integers
 * @param N Number of rows
 * @param M Number of columns
 * @return Pointer to allocated array, initialized to zero
 */
int **i_2d_allocate(const long unsigned int N, const long unsigned int M);

/**
 * @brief Free 2D array of integers
 * @param A Pointer to array to free
 */
void free_i_2d_allocate(int **A);

/**
 * @brief Allocate 3D array of integers
 * @param N First dimension size
 * @param M Second dimension size
 * @param L Third dimension size
 * @return Pointer to allocated array, initialized to zero
 */
int ***i_3d_allocate(const long unsigned int N, const long unsigned int M, const long unsigned int L);

/**
 * @brief Free 3D array of integers
 * @param A Pointer to array to free
 */
void free_i_3d_allocate(int ***A);

/**
 * @brief Allocate 1D array of doubles
 * @param N Size of array to allocate
 * @return Pointer to allocated array, initialized to zero
 */
double *d_1d_allocate(const long unsigned int N);

/**
 * @brief Free 1D array of doubles
 * @param A Pointer to array to free
 */
void free_d_1d_allocate(double *A);

/**
 * @brief Allocate 2D array of doubles
 * @param N Number of rows
 * @param M Number of columns
 * @return Pointer to allocated array, initialized to zero
 */
double **d_2d_allocate(const long unsigned int N, const long unsigned int M);

/**
 * @brief Free 2D array of doubles
 * @param A Pointer to array to free
 */
void free_d_2d_allocate(double **A);

/**
 * @brief Allocate 1D array of complex doubles
 * @param N Size of array to allocate
 * @return Pointer to allocated array, initialized to zero
 */
complex double *cd_1d_allocate(const long unsigned int N);

/**
 * @brief Free 1D array of complex doubles
 * @param A Pointer to array to free
 */
void free_cd_1d_allocate(double complex*A);

/**
 * @brief Allocate 2D array of complex doubles
 * @param N Number of rows
 * @param M Number of columns
 * @return Pointer to allocated array, initialized to zero
 */
complex double **cd_2d_allocate(const long unsigned int N, const long unsigned int M);

/**
 * @brief Free 2D array of complex doubles
 * @param A Pointer to array to free
 */
void free_cd_2d_allocate(double complex**A);

/**
 * @brief Allocate 3D array of complex doubles
 * @param N First dimension size
 * @param M Second dimension size
 * @param L Third dimension size
 * @return Pointer to allocated array, initialized to zero
 */
double complex***cd_3d_allocate(const long unsigned int N, const long unsigned int M, const long unsigned int L);

/**
 * @brief Free 3D array of complex doubles
 * @param A Pointer to array to free
 */
void free_cd_3d_allocate(double complex***A);

#endif //MVMC_SETMEMORY_H
