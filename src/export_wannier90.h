/**
 * @file export_wannier90.h
 * @brief Header for exporting model data in Wannier90 format
 *
 * @details
 * Declares functions to export lattice geometry and interaction parameters
 * in a format compatible with Wannier90. These functions are only available
 * when compiled with the _HWAVE flag.
 *
 * @copyright
 * HPhi-mVMC-StdFace - Common input generator
 * Copyright (C) 2015 The University of Tokyo
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */
#ifndef EXPORT_WANNIER90_INCLUDED
#define EXPORT_WANNIER90_INCLUDED

//#include "StdFace_vals.h"

#if defined(_HWAVE)

/**
 * @brief Export geometry data in Wannier90 format
 *
 * @details
 * Writes lattice vectors and orbital positions to a file in the format
 * required by Wannier90. The output includes primitive lattice vectors,
 * number of sites per unit cell, and orbital positions in fractional
 * coordinates.
 *
 * @param[in] StdI Standard input parameters containing lattice geometry
 *                  information including primitive vectors, number of sites
 *                  per unit cell, and orbital positions
 */
void ExportGeometry(struct StdIntList *StdI);

/**
 * @brief Export interaction parameters in Wannier90 format
 *
 * @details
 * Writes interaction parameters between orbitals to files in the format
 * required by Wannier90. Exports transfer, Coulomb (intra and inter),
 * Hund, exchange, pair-lift, and pair-hopping coefficients.
 *
 * @param[in] StdI Standard input parameters containing interaction terms
 *                  between orbitals, including hopping and correlation
 *                  parameters
 */
void ExportInteraction(struct StdIntList *StdI);

#endif

#endif /* EXPORT_WANNIER90_INCLUDED */
