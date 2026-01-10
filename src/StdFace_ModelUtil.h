/*
HPhi-mVMC-StdFace - Common input generator
Copyright (C) 2015 The University of Tokyo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @file StdFace_ModelUtil.h
 * @brief Utility functions for constructing models in standard mode
 * @author Mitsuaki Kawamura (The University of Tokyo)
 * 
 * This file contains declarations for utility functions used in constructing
 * various lattice models and handling interactions in standard mode.
 */

#include <complex.h>
#include <stdio.h>

/**
 * @brief MPI Abortation wrapper
 * @param[in] errorcode Error code to exit with
 */
void StdFace_exit(int errorcode);

/**
 * @brief Add interaction term to the list
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] intr0 Interaction strength
 * @param[in] site1 First site index
 * @param[in] spin1 First spin index
 * @param[in] site2 Second site index
 * @param[in] spin2 Second spin index
 * @param[in] site3 Third site index
 * @param[in] spin3 Third spin index
 * @param[in] site4 Fourth site index
 * @param[in] spin4 Fourth spin index
 */
void StdFace_intr(struct StdIntList *StdI, double complex intr0,
  int site1, int spin1, int site2, int spin2,
  int site3, int spin3, int site4, int spin4);

/**
 * @brief Add hopping term between sites
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] trans0 Hopping amplitude
 * @param[in] isite Source site index
 * @param[in] jsite Target site index
 * @param[in] dR Distance vector between sites
 */
void StdFace_Hopping(struct StdIntList *StdI, double complex trans0, int isite, int jsite, double *dR);

/**
 * @brief Add transfer term between sites with spin
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] trans0 Transfer amplitude
 * @param[in] isite Source site index
 * @param[in] ispin Source spin index
 * @param[in] jsite Target site index
 * @param[in] jspin Target spin index
 */
void StdFace_trans(struct StdIntList *StdI,double complex trans0,int isite,int ispin,int jsite,int jspin);

/**
 * @brief Add local Hubbard terms
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] mu0 Chemical potential
 * @param[in] h0 Magnetic field
 * @param[in] Gamma0 Gamma field strength
 * @param[in] Gamma0_y Y-component of Gamma field
 * @param[in] U0 On-site Coulomb interaction
 * @param[in] isite Site index
 */
void StdFace_HubbardLocal(struct StdIntList *StdI, double mu0, double h0,
  double Gamma0, double Gamma0_y, double U0, int isite);

/**
 * @brief Add magnetic field terms
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] S2 Twice the spin value
 * @param[in] h Magnetic field strength
 * @param[in] Gamma Gamma field strength
 * @param[in] Gamma_y Y-component of Gamma field
 * @param[in] isite Site index
 */
void StdFace_MagField(struct StdIntList *StdI, int S2, double h, double Gamma, double Gamma_y, int isite);

/**
 * @brief Add Coulomb interaction between sites
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] V Coulomb interaction strength
 * @param[in] isite First site index
 * @param[in] jsite Second site index
 */
void StdFace_Coulomb(struct StdIntList *StdI, double V, int isite, int jsite);

/**
 * @brief Add general exchange interaction
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] J Exchange coupling matrix
 * @param[in] Si2 Twice the spin value for site i
 * @param[in] Sj2 Twice the spin value for site j
 * @param[in] isite First site index
 * @param[in] jsite Second site index
 */
void StdFace_GeneralJ(struct StdIntList *StdI, double J[3][3],
  int Si2, int Sj2, int isite, int jsite);

/**
 * @brief Print double value with name
 * @param[in] valname Name of the value
 * @param[in] val Pointer to value
 * @param[in] val0 Default value
 */
void StdFace_PrintVal_d(char* valname, double *val, double val0);

/**
 * @brief Print two double values with name
 * @param[in] valname Name of the values
 * @param[in] val Pointer to values
 * @param[in] val0 First default value
 * @param[in] val1 Second default value
 */
void StdFace_PrintVal_dd(char* valname, double *val, double val0, double val1);

/**
 * @brief Print complex value with name
 * @param[in] valname Name of the value
 * @param[in] val Pointer to value
 * @param[in] val0 Default value
 */
void StdFace_PrintVal_c(char* valname, double complex *val, double complex val0);

/**
 * @brief Print integer value with name
 * @param[in] valname Name of the value
 * @param[in] val Pointer to value
 * @param[in] val0 Default value
 */
void StdFace_PrintVal_i(char* valname, int *val, int val0);

/**
 * @brief Print warning for unused double parameter
 * @param[in] valname Name of unused parameter
 * @param[in] val Value of unused parameter
 */
void StdFace_NotUsed_d(char* valname, double val);

/**
 * @brief Print warning for unused integer parameter
 * @param[in] valname Name of unused parameter
 * @param[in] val Value of unused parameter
 */
void StdFace_NotUsed_i(char* valname, int val);

/**
 * @brief Print warning for unused complex parameter
 * @param[in] valname Name of unused parameter
 * @param[in] val Value of unused parameter
 */
void StdFace_NotUsed_c(char* valname, double complex val);

/**
 * @brief Print warning for unused exchange parameter
 * @param[in] valname Name of unused parameter
 * @param[in] JAll Scalar exchange value
 * @param[in] J Exchange matrix
 */
void StdFace_NotUsed_J(char* valname, double JAll, double J[3][3]);

/**
 * @brief Check if required integer parameter is set
 * @param[in] valname Name of parameter
 * @param[in] val Value to check
 */
void StdFace_RequiredVal_i(char* valname, int val);

/**
 * @brief Initialize spin nearest neighbor interactions
 * @param[out] J Output exchange matrix
 * @param[in] JAll Scalar exchange value
 * @param[in] J0 Input exchange matrix
 * @param[in] J0All Input scalar exchange value
 * @param[in] J0name Name of exchange parameter
 */
void StdFace_InputSpinNN(double J[3][3], double JAll, double J0[3][3], double J0All, char *J0name);

/**
 * @brief Initialize spin interactions
 * @param[out] Jp Output exchange matrix
 * @param[in] JpAll Scalar exchange value
 * @param[in] Jpname Name of exchange parameter
 */
void StdFace_InputSpin(double Jp[3][3], double JpAll, char *Jpname);

/**
 * @brief Initialize Coulomb interaction
 * @param[in] V Input Coulomb strength
 * @param[out] V0 Output Coulomb strength
 * @param[in] V0name Name of Coulomb parameter
 */
void StdFace_InputCoulombV(double V, double *V0, char *V0name);

/**
 * @brief Initialize hopping parameter
 * @param[in] t Input hopping
 * @param[out] t0 Output hopping
 * @param[in] t0name Name of hopping parameter
 */
void StdFace_InputHopp(double complex t, double complex *t0, char *t0name);

/**
 * @brief Initialize lattice sites
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] fp File pointer for output
 * @param[in] dim Lattice dimension
 */
void StdFace_InitSite(struct StdIntList *StdI, FILE *fp, int dim);

/**
 * @brief Set site labels and connections
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] fp File pointer for output
 * @param[in] iW Width index
 * @param[in] iL Length index
 * @param[in] diW Width displacement
 * @param[in] diL Length displacement
 * @param[in] isiteUC First unit cell site
 * @param[in] jsiteUC Second unit cell site
 * @param[out] isite Output first site index
 * @param[out] jsite Output second site index
 * @param[in] connect Connection type
 * @param[out] Cphase Complex phase factor
 * @param[out] dR Distance vector
 */
void StdFace_SetLabel(struct StdIntList *StdI, FILE *fp,
  int iW, int iL, int diW, int diL, int isiteUC, int jsiteUC,
  int *isite, int *jsite, int connect, double complex *Cphase, double *dR);

/**
 * @brief Print geometry information
 * @param[in] StdI Structure containing model parameters
 */
void StdFace_PrintGeometry(struct StdIntList *StdI);

/**
 * @brief Allocate memory for interactions
 * @param[in,out] StdI Structure containing model parameters
 * @param[in] ntransMax Maximum number of transfer terms
 * @param[in] nintrMax Maximum number of interaction terms
 */
void StdFace_MallocInteractions(struct StdIntList *StdI, int ntransMax, int nintrMax);

/**
 * @brief Find site indices in lattice
 * @param[in] StdI Structure containing model parameters
 * @param[in] iW Width index
 * @param[in] iL Length index
 * @param[in] iH Height index
 * @param[in] diW Width displacement
 * @param[in] diL Length displacement
 * @param[in] diH Height displacement
 * @param[in] isiteUC First unit cell site
 * @param[in] jsiteUC Second unit cell site
 * @param[out] isite Output first site index
 * @param[out] jsite Output second site index
 * @param[out] Cphase Complex phase factor
 * @param[out] dR Distance vector
 */
void StdFace_FindSite(struct StdIntList *StdI,
  int iW, int iL, int iH, int diW, int diL, int diH,
  int isiteUC, int jsiteUC,
  int *isite, int *jsite, double complex *Cphase, double *dR);

/**
 * @brief Print XSF format output
 * @param[in] StdI Structure containing model parameters
 */
void StdFace_PrintXSF(struct StdIntList *StdI);

/**
 * @brief Initialize tetragonal lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Tetragonal(struct StdIntList *StdI);

/**
 * @brief Initialize chain lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Chain(struct StdIntList *StdI);

/**
 * @brief Initialize ladder lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Ladder(struct StdIntList *StdI);

/**
 * @brief Initialize triangular lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Triangular(struct StdIntList *StdI);

/**
 * @brief Initialize honeycomb lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Honeycomb(struct StdIntList *StdI);

/**
 * @brief Initialize kagome lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Kagome(struct StdIntList *StdI);

/**
 * @brief Initialize orthorhombic lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Orthorhombic(struct StdIntList *StdI);

/**
 * @brief Initialize face-centered orthorhombic lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_FCOrtho(struct StdIntList *StdI);

/**
 * @brief Initialize pyrochlore lattice
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Pyrochlore(struct StdIntList *StdI);

/**
 * @brief Initialize from Wannier90 input
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Wannier90(struct StdIntList *StdI);

#if defined(_HPhi)
/**
 * @brief Initialize chain lattice with boost
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Chain_Boost(struct StdIntList *StdI);

/**
 * @brief Initialize ladder lattice with boost
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Ladder_Boost(struct StdIntList *StdI);

/**
 * @brief Initialize honeycomb lattice with boost
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Honeycomb_Boost(struct StdIntList *StdI);

/**
 * @brief Initialize kagome lattice with boost
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Kagome_Boost(struct StdIntList *StdI);
#elif defined(_mVMC)
/**
 * @brief Generate orbital information
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_generate_orb(struct StdIntList *StdI);

/**
 * @brief Handle projections
 * @param[in,out] StdI Structure containing model parameters
 */
void StdFace_Proj(struct StdIntList *StdI);

/**
 * @brief Print Jastrow factors
 * @param[in] StdI Structure containing model parameters
 */
void PrintJastrow(struct StdIntList *StdI);
#endif
