/**
 * @file Ladder.c
 * @brief Standard mode for the Ladder lattice
 * @author Mitsuaki Kawamura (The University of Tokyo)
 *
 * @details
 * This file contains functions to set up Hamiltonians for ladder lattice models.
 * It supports Heisenberg, Hubbard and Kondo models with various interactions.
 *
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

#include "StdFace_vals.h"
#include "StdFace_ModelUtil.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

/**
 * @brief Setup a Hamiltonian for the generalized Heisenberg model on a ladder lattice
 *
 * @details
 * This function sets up the Hamiltonian parameters for a ladder lattice model.
 * It supports:
 * - Heisenberg model with spin interactions
 * - Hubbard model with electron hopping and interactions
 * - Kondo model combining localized spins and itinerant electrons
 *
 * The ladder geometry consists of:
 * - Vertical rungs between the two chains
 * - Nearest and next-nearest neighbor interactions along the chains
 * - Diagonal interactions between the chains
 *
 * The function:
 * 1. Sets up the lattice geometry and parameters
 * 2. Validates input parameters
 * 3. Allocates arrays for interactions
 * 4. Sets up all the interaction terms in the Hamiltonian
 *
 * @param[in,out] StdI Pointer to the structure containing model parameters
 */
void StdFace_Ladder(
  struct StdIntList *StdI
)
{
  FILE *fp = NULL;
  int isite, jsite, ntransMax, nintrMax;
  int iL, isiteUC;
  double complex Cphase;
  double dR[3];

  /* Open file for lattice plot if enabled */
#ifdef _HWAVE
  if (StdI->lattice_gp == 1)
#endif
  fp = fopen("lattice.gp", "w");

  /*
   * 1. Set lattice size and shape parameters
   */
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  StdFace_PrintVal_d("Wlength", &StdI->length[0], StdI->a);
  StdFace_PrintVal_d("Llength", &StdI->length[1], StdI->a);
  StdFace_PrintVal_d("Wx", &StdI->direct[0][0], StdI->length[0]);
  StdFace_PrintVal_d("Wy", &StdI->direct[0][1], 0.0);
  StdFace_PrintVal_d("Lx", &StdI->direct[1][0], 0.0);
  StdFace_PrintVal_d("Ly", &StdI->direct[1][1], StdI->length[1]);

  /* Required parameters */
  StdFace_RequiredVal_i("L", StdI->L);
  StdFace_RequiredVal_i("W", StdI->W);
  
  /* Unused parameters */
  StdFace_NotUsed_i("a0W", StdI->box[0][0]);
  StdFace_NotUsed_i("a0L", StdI->box[0][1]);
  StdFace_NotUsed_i("a1W", StdI->box[1][0]);
  StdFace_NotUsed_i("a1L", StdI->box[1][1]);

  /* Phase factors */
  StdFace_PrintVal_d("phase0", &StdI->phase[0], 0.0);
  StdFace_NotUsed_d("phase1", StdI->phase[1]);
  StdI->phase[1] = StdI->phase[0];
  StdI->phase[0] = 0.0;

  /* Set unit cell */
  StdI->NsiteUC = StdI->W;
  StdI->W = 1;
  StdI->direct[0][0] = (double)StdI->NsiteUC;
  StdFace_InitSite(StdI, fp, 2);
  
  /* Initialize site positions */
  for (isite = 0; isite < StdI->NsiteUC; isite++){
    StdI->tau[isite][0] = (double)isite / (double)StdI->NsiteUC;
    StdI->tau[isite][1] = 0.0; 
    StdI->tau[isite][2] = 0.0;
  }

  /*
   * 2. Set Hamiltonian parameters
   */
  fprintf(stdout, "\n  @ Hamiltonian \n\n");

  /* Unused coupling parameters */
  StdFace_NotUsed_J("J", StdI->JAll, StdI->J);
  StdFace_NotUsed_J("J'", StdI->JpAll, StdI->Jp);
  StdFace_NotUsed_c("t", StdI->t);
  StdFace_NotUsed_c("t'", StdI->tp);
  StdFace_NotUsed_d("V", StdI->V);
  StdFace_NotUsed_d("V'", StdI->Vp);
  StdFace_NotUsed_d("K", StdI->K);

  /* Magnetic field parameters */
  StdFace_PrintVal_d("h", &StdI->h, 0.0);
  StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
  StdFace_PrintVal_d("Gamma_y", &StdI->Gamma_y, 0.0);

  /* Model specific parameters */
  if (strcmp(StdI->model, "spin") == 0 ) {
    /* Heisenberg model parameters */
    StdFace_PrintVal_i("2S", &StdI->S2, 1);
    StdFace_PrintVal_d("D", &StdI->D[2][2], 0.0);
    StdFace_InputSpin(StdI->J0, StdI->J0All, "J0");
    StdFace_InputSpin(StdI->J1, StdI->J1All, "J1");
    StdFace_InputSpin(StdI->J2, StdI->J2All, "J2");
    StdFace_InputSpin(StdI->J1p, StdI->J1pAll, "J1'");
    StdFace_InputSpin(StdI->J2p, StdI->J2pAll, "J2'");

    /* Unused electronic parameters */
    StdFace_NotUsed_d("mu", StdI->mu);
    StdFace_NotUsed_d("U", StdI->U);
    StdFace_NotUsed_c("t0", StdI->t0);
    StdFace_NotUsed_c("t1", StdI->t1);
    StdFace_NotUsed_c("t2", StdI->t2);
    StdFace_NotUsed_c("t1'", StdI->t1p);
    StdFace_NotUsed_c("t2'", StdI->t2p);
    StdFace_NotUsed_d("V0", StdI->V0);
    StdFace_NotUsed_d("V1", StdI->V1);
    StdFace_NotUsed_d("V2", StdI->V2);
    StdFace_NotUsed_d("V1'", StdI->V1p);
    StdFace_NotUsed_d("V2'", StdI->V2p);

  } else {
    /* Electronic model parameters */
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_InputHopp(StdI->t, &StdI->t0, "t0");
    StdFace_InputHopp(StdI->t, &StdI->t1, "t1");
    StdFace_InputHopp(StdI->t, &StdI->t2, "t2");
    StdFace_InputHopp(StdI->t, &StdI->t1p, "t1'");
    StdFace_InputHopp(StdI->t, &StdI->t2p, "t2'");
    StdFace_InputCoulombV(StdI->V, &StdI->V0, "V0");
    StdFace_InputCoulombV(StdI->V, &StdI->V1, "V1");
    StdFace_InputCoulombV(StdI->V, &StdI->V2, "V2");
    StdFace_InputCoulombV(StdI->V, &StdI->V1p, "V1'");
    StdFace_InputCoulombV(StdI->V, &StdI->V2p, "V2'");

    /* Unused spin parameters */
    StdFace_NotUsed_J("J0", StdI->J0All, StdI->J0);
    StdFace_NotUsed_J("J1", StdI->J1All, StdI->J1);
    StdFace_NotUsed_J("J2", StdI->J2All, StdI->J2);
    StdFace_NotUsed_J("J1p", StdI->J1pAll, StdI->J1p);
    StdFace_NotUsed_J("J2p", StdI->J2pAll, StdI->J2p);
    StdFace_NotUsed_d("D", StdI->D[2][2]);

    if (strcmp(StdI->model, "hubbard") == 0 ) {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_J("J", StdI->JAll, StdI->J);
    }
    else {
      StdFace_PrintVal_i("2S", &StdI->S2, 1);
      StdFace_InputSpin(StdI->J, StdI->JAll, "J");
    }
  }

  fprintf(stdout, "\n  @ Numerical conditions\n\n");

  /*
   * 3. Set local spin flags and number of sites
   */
  StdI->nsite = StdI->L * StdI->NsiteUC;
  if (strcmp(StdI->model, "kondo") == 0 ) StdI->nsite *= 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  
  /* Set local spin flags based on model type */
  if (strcmp(StdI->model, "spin") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = StdI->S2;
  else if (strcmp(StdI->model, "hubbard") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++)StdI->locspinflag[isite] = 0;
  else if (strcmp(StdI->model, "kondo") == 0 )
    for (isite = 0; isite < StdI->nsite / 2; isite++) {
      StdI->locspinflag[isite] = StdI->S2;
      StdI->locspinflag[isite + StdI->nsite / 2] = 0;
    }

  /*
   * 4. Calculate maximum number of interactions and allocate arrays
   */
  if (strcmp(StdI->model, "spin") == 0 ) {
    /* For spin model */
    ntransMax = StdI->L * StdI->NsiteUC * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    nintrMax = StdI->L * StdI->NsiteUC * (1/*D*/ + 1/*J1*/ + 1/*J1'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1)
      + StdI->L * (StdI->NsiteUC - 1) * (1/*J0*/ + 1/*J2*/ + 1/*J2'*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  } else {
    /* For electronic models */
    ntransMax = StdI->L*StdI->NsiteUC * 2/*spin*/ * (2/*mu+h+Gamma*/ + 2/*t1*/ + 2/*t1'*/)
      + StdI->L*(StdI->NsiteUC - 1) * 2/*spin*/ * (2/*t0*/ + 2/*t2*/ + 2/*t2'*/);
    nintrMax = StdI->L*StdI->NsiteUC * 1/*U*/
      + StdI->L*StdI->NsiteUC * 4 * (1/*V1*/ + 1/*V1'*/)
      + StdI->L*(StdI->NsiteUC - 1) * 4 * (1/*V0*/ + 1/*V2*/ + 1/*V2'*/);

    if (strcmp(StdI->model, "kondo") == 0) {
      ntransMax += StdI->L * StdI->NsiteUC * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
      nintrMax += StdI->nsite / 2 * (3 * 1 + 1) * (3 * StdI->S2 + 1);
    }
  }

  /* Allocate arrays */
  StdFace_MallocInteractions(StdI, ntransMax, nintrMax);

  /*
   * 5. Set all interactions
   */
  for (iL = 0; iL < StdI->L; iL++) {
    for (isiteUC = 0; isiteUC < StdI->NsiteUC; isiteUC++) {

      isite = isiteUC + iL * StdI->NsiteUC;
      if (strcmp(StdI->model, "kondo") == 0 ) isite += StdI->L * StdI->NsiteUC;

      /* Local terms */
      if (strcmp(StdI->model, "spin") == 0 ) {
        StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, -StdI->Gamma_y, isite);
        StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite, isite);
      } else {
        StdFace_HubbardLocal(StdI, StdI->mu, -StdI->h, -StdI->Gamma, -StdI->Gamma_y, StdI->U, isite);
        if (strcmp(StdI->model, "kondo") == 0 ) {
          jsite = isiteUC + iL * StdI->NsiteUC;
          StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
          StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, -StdI->Gamma_y, jsite);
        }
      }

      /* Nearest neighbor along the ladder */
      StdFace_SetLabel(StdI, fp, 0, iL, 0, 1, isiteUC, isiteUC, &isite, &jsite, 1, &Cphase, dR);
      
      if (strcmp(StdI->model, "spin") == 0 ) {
        StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
      } else {
        StdFace_Hopping(StdI, Cphase * StdI->t1, isite, jsite, dR);
        StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
      }

      /* Second nearest neighbor along the ladder */
      StdFace_SetLabel(StdI, fp, 0, iL, 0, 2, isiteUC, isiteUC, &isite, &jsite, 2, &Cphase, dR);
      
      if (strcmp(StdI->model, "spin") == 0 ) {
        StdFace_GeneralJ(StdI, StdI->J1p, StdI->S2, StdI->S2, isite, jsite);
      } else {
        StdFace_Hopping(StdI, Cphase * StdI->t1p, isite, jsite, dR);
        StdFace_Coulomb(StdI, StdI->V1p, isite, jsite);
      }

      /* Interactions across rungs */
      if (isiteUC < StdI->NsiteUC - 1) {
        /* Vertical */
        StdFace_SetLabel(StdI, fp, 0, iL, 0, 0, isiteUC, isiteUC + 1, &isite, &jsite, 1, &Cphase, dR);
        
        if (strcmp(StdI->model, "spin") == 0 ) {
          StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
        } else {
          StdFace_Hopping(StdI, Cphase * StdI->t0, isite, jsite, dR);
          StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
        }

        /* Diagonal 1 */
        StdFace_SetLabel(StdI, fp, 0, iL, 0, 1, isiteUC, isiteUC + 1, &isite, &jsite, 1, &Cphase, dR);
        
        if (strcmp(StdI->model, "spin") == 0 ) {
          StdFace_GeneralJ(StdI, StdI->J2, StdI->S2, StdI->S2, isite, jsite);
        } else {
          StdFace_Hopping(StdI, Cphase * StdI->t2, isite, jsite, dR);
          StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
        }

        /* Diagonal 2 */
        StdFace_SetLabel(StdI, fp, 0, iL, 0, -1, isiteUC, isiteUC + 1, &isite, &jsite, 1, &Cphase, dR);
        
        if (strcmp(StdI->model, "spin") == 0 ) {
          StdFace_GeneralJ(StdI, StdI->J2p, StdI->S2, StdI->S2, isite, jsite);
        } else {
          StdFace_Hopping(StdI, Cphase * StdI->t2p, isite, jsite, dR);
          StdFace_Coulomb(StdI, StdI->V2p, isite, jsite);
        }

      }
    }
  }

  /* Close lattice plot file */
#ifdef _HWAVE
  if (StdI->lattice_gp == 1) {
#endif
  fprintf(fp, "plot \'-\' w d lc 7\n0.0 0.0\nend\npause -1\n");
  fclose(fp);
#ifdef _HWAVE
  }
#endif

  /* Print final geometry information */
  StdFace_PrintGeometry(StdI);
}

#if defined(_HPhi)
/**
 * @brief Setup a Hamiltonian for the generalized Heisenberg model on a ladder lattice with boost
 *
 * @details
 * This function sets up a boosted version of the ladder Hamiltonian for HPhi.
 * It is specialized for S=1/2 Heisenberg models and requires W=2, even L>=4.
 * The boost method writes interaction parameters and 6-spin pair lists to
 * a "boost.def" file used by the HPhi solver.
 *
 * @param[in,out] StdI Pointer to the structure containing model parameters
 */
void StdFace_Ladder_Boost(struct StdIntList *StdI)
{
  int isite, ipivot;
  int kintr;
  FILE *fp;

  /* Validate and set parameters */
  StdI->W = StdI->NsiteUC;
  StdI->NsiteUC = 1;

  /* Open boost definition file */
  fp = fopen("boost.def", "w");

  /* Write magnetic field parameters */
  fprintf(fp, "# Magnetic field\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    -0.5 * StdI->Gamma, -0.5 * StdI->Gamma_y, -0.5 * StdI->h);

  /* Write interaction parameters */
  fprintf(fp, "%d  # Number of type of J\n", 5);

  /* J1 - Vertical interactions */
  fprintf(fp, "# J 1 (inter chain, vertical)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][0], 0.25 * StdI->J0[0][1], 0.25 * StdI->J0[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n", 
    0.25 * StdI->J0[0][1], 0.25 * StdI->J0[1][1], 0.25 * StdI->J0[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J0[0][2], 0.25 * StdI->J0[1][2], 0.25 * StdI->J0[2][2]);

  /* J2 - Nearest neighbor along chain */
  fprintf(fp, "# J 2 (Nearest neighbor, along chain)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][0], 0.25 * StdI->J1[0][1], 0.25 * StdI->J1[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1[0][1], 0.25 * StdI->J1[1][1], 0.25 * StdI->J1[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n", 
    0.25 * StdI->J1[0][2], 0.25 * StdI->J1[1][2], 0.25 * StdI->J1[2][2]);

  /* J3 - Second nearest neighbor along chain */
  fprintf(fp, "# J 3 (Second nearest neighbor, along chain)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1p[0][0], 0.25 * StdI->J1p[0][1], 0.25 * StdI->J1p[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1p[0][1], 0.25 * StdI->J1p[1][1], 0.25 * StdI->J1p[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J1p[0][2], 0.25 * StdI->J1p[1][2], 0.25 * StdI->J1p[2][2]);

  /* J4 - Diagonal 1 interactions */
  fprintf(fp, "# J 4 (inter chain, diagonal1)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][0], 0.25 * StdI->J2[0][1], 0.25 * StdI->J2[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][1], 0.25 * StdI->J2[1][1], 0.25 * StdI->J2[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2[0][2], 0.25 * StdI->J2[1][2], 0.25 * StdI->J2[2][2]);

  /* J5 - Diagonal 2 interactions */
  fprintf(fp, "# J 5 (inter chain, diagonal2)\n");
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2p[0][0], 0.25 * StdI->J2p[0][1], 0.25 * StdI->J2p[0][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2p[0][1], 0.25 * StdI->J2p[1][1], 0.25 * StdI->J2p[1][2]);
  fprintf(fp, "%25.15e %25.15e %25.15e\n",
    0.25 * StdI->J2p[0][2], 0.25 * StdI->J2p[1][2], 0.25 * StdI->J2p[2][2]);

  /* Validate parameters */
  if (StdI->S2 != 1) {
    fprintf(stdout, "\n ERROR! S2 must be 1 in Boost. \n\n");
    StdFace_exit(-1);
  }
  StdI->ishift_nspin = 2;
  if (StdI->W != 2) {
    fprintf(stdout, "\n ERROR! W != 2 \n\n");
    StdFace_exit(-1);
  }
  if (StdI->L % 2 != 0) {
    fprintf(stdout, "\n ERROR! L %% 2 != 0 \n\n");
    StdFace_exit(-1);
  }
  if (StdI->L < 4) {
    fprintf(stdout, "\n ERROR! L < 4 \n\n");
    StdFace_exit(-1);
  }

  /* Set dimensions */
  StdI->W = StdI->L;
  StdI->L = 2;
  StdI->num_pivot = StdI->W / 2;

  /* Write topology information */
  fprintf(fp, "# W0  R0  StdI->num_pivot  StdI->ishift_nspin\n");
  fprintf(fp, "%d %d %d %d\n", StdI->W, StdI->L, StdI->num_pivot, StdI->ishift_nspin);

  /* Allocate and initialize 6-spin star list */
  StdI->list_6spin_star = (int **)malloc(sizeof(int*) * StdI->num_pivot);
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_star[ipivot] = (int *)malloc(sizeof(int) * 7);
  }

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_star[ipivot][0] = 7; // num of J
    StdI->list_6spin_star[ipivot][1] = 1;
    StdI->list_6spin_star[ipivot][2] = 1;
    StdI->list_6spin_star[ipivot][3] = 1;
    StdI->list_6spin_star[ipivot][4] = 1;
    StdI->list_6spin_star[ipivot][5] = 1;
    StdI->list_6spin_star[ipivot][6] = 1; // flag
  }

  fprintf(fp, "# StdI->list_6spin_star\n");
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    fprintf(fp, "# pivot %d\n", ipivot);
    for (isite = 0; isite < 7; isite++) {
      fprintf(fp, "%d ", StdI->list_6spin_star[ipivot][isite]);
    }
    fprintf(fp, "\n");
  }

  StdI->list_6spin_pair = (int ***)malloc(sizeof(int**) * StdI->num_pivot);
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_pair[ipivot] = (int **)malloc(sizeof(int*) * 7);
    for (isite = 0; isite < 7; isite++) {
      StdI->list_6spin_pair[ipivot][isite] = (int *)malloc(sizeof(int) * StdI->list_6spin_star[ipivot][0]);
    }
  }

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    StdI->list_6spin_pair[ipivot][0][0] = 0;
    StdI->list_6spin_pair[ipivot][1][0] = 1;
    StdI->list_6spin_pair[ipivot][2][0] = 2;
    StdI->list_6spin_pair[ipivot][3][0] = 3;
    StdI->list_6spin_pair[ipivot][4][0] = 4;
    StdI->list_6spin_pair[ipivot][5][0] = 5;
    StdI->list_6spin_pair[ipivot][6][0] = 1; // type of J
    StdI->list_6spin_pair[ipivot][0][1] = 0;
    StdI->list_6spin_pair[ipivot][1][1] = 2;
    StdI->list_6spin_pair[ipivot][2][1] = 1;
    StdI->list_6spin_pair[ipivot][3][1] = 3;
    StdI->list_6spin_pair[ipivot][4][1] = 4;
    StdI->list_6spin_pair[ipivot][5][1] = 5;
    StdI->list_6spin_pair[ipivot][6][1] = 2; // type of J
    StdI->list_6spin_pair[ipivot][0][2] = 1;
    StdI->list_6spin_pair[ipivot][1][2] = 3;
    StdI->list_6spin_pair[ipivot][2][2] = 0;
    StdI->list_6spin_pair[ipivot][3][2] = 2;
    StdI->list_6spin_pair[ipivot][4][2] = 4;
    StdI->list_6spin_pair[ipivot][5][2] = 5;
    StdI->list_6spin_pair[ipivot][6][2] = 2; // type of J
    StdI->list_6spin_pair[ipivot][0][3] = 0;
    StdI->list_6spin_pair[ipivot][1][3] = 4;
    StdI->list_6spin_pair[ipivot][2][3] = 1;
    StdI->list_6spin_pair[ipivot][3][3] = 2;
    StdI->list_6spin_pair[ipivot][4][3] = 3;
    StdI->list_6spin_pair[ipivot][5][3] = 5;
    StdI->list_6spin_pair[ipivot][6][3] = 3; // type of J
    StdI->list_6spin_pair[ipivot][0][4] = 1;
    StdI->list_6spin_pair[ipivot][1][4] = 5;
    StdI->list_6spin_pair[ipivot][2][4] = 0;
    StdI->list_6spin_pair[ipivot][3][4] = 2;
    StdI->list_6spin_pair[ipivot][4][4] = 3;
    StdI->list_6spin_pair[ipivot][5][4] = 4;
    StdI->list_6spin_pair[ipivot][6][4] = 3; // type of J
    StdI->list_6spin_pair[ipivot][0][5] = 0;
    StdI->list_6spin_pair[ipivot][1][5] = 3;
    StdI->list_6spin_pair[ipivot][2][5] = 1;
    StdI->list_6spin_pair[ipivot][3][5] = 2;
    StdI->list_6spin_pair[ipivot][4][5] = 4;
    StdI->list_6spin_pair[ipivot][5][5] = 5;
    StdI->list_6spin_pair[ipivot][6][5] = 4; // type of J
    StdI->list_6spin_pair[ipivot][0][6] = 1;
    StdI->list_6spin_pair[ipivot][1][6] = 2;
    StdI->list_6spin_pair[ipivot][2][6] = 0;
    StdI->list_6spin_pair[ipivot][3][6] = 3;
    StdI->list_6spin_pair[ipivot][4][6] = 4;
    StdI->list_6spin_pair[ipivot][5][6] = 5;
    StdI->list_6spin_pair[ipivot][6][6] = 5; // type of J
  }

  fprintf(fp, "# StdI->list_6spin_pair\n");
  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    fprintf(fp, "# pivot %d\n", ipivot);
    for (kintr = 0; kintr < StdI->list_6spin_star[ipivot][0]; kintr++) {
      for (isite = 0; isite < 7; isite++) {
        fprintf(fp, "%d ", StdI->list_6spin_pair[ipivot][isite][kintr]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    free(StdI->list_6spin_star[ipivot]);
  }
  free(StdI->list_6spin_star);

  for (ipivot = 0; ipivot < StdI->num_pivot; ipivot++) {
    for (isite = 0; isite < 7; isite++) {
      free(StdI->list_6spin_pair[ipivot][isite]);
    }
    free(StdI->list_6spin_pair[ipivot]);
  }
  free(StdI->list_6spin_pair);
}
#endif
