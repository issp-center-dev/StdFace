#ifndef EXPORT_WANNIER90_INCLUDED
#define EXPORT_WANNIER90_INCLUDED

//#include "StdFace_vals.h"

#if defined(_HWAVE)

/**
 * Export geometry data in Wannier90 format
 *
 * Parameters
 * ----------
 * StdI : struct StdIntList*
 *     Standard input parameters containing lattice geometry information
 *     including primitive vectors, number of sites per unit cell,
 *     and orbital positions
 *
 * Notes
 * -----
 * Writes lattice vectors and orbital positions to a file in the format
 * required by Wannier90
 */
void ExportGeometry(struct StdIntList *StdI);

/**
 * Export interaction parameters in Wannier90 format
 *
 * Parameters
 * ----------
 * StdI : struct StdIntList*
 *     Standard input parameters containing interaction terms between
 *     orbitals, including hopping and correlation parameters
 *
 * Notes
 * -----
 * Writes interaction parameters between orbitals to a file in the
 * format required by Wannier90
 */
void ExportInteraction(struct StdIntList *StdI);

#endif

#endif /* EXPORT_WANNIER90_INCLUDED */
