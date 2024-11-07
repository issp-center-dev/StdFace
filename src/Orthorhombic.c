/**
 * @file Orthorhombic.c
 * @brief Standard mode for the orthorhombic lattice
 * @copyright Copyright (C) 2015 The University of Tokyo
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
 * @brief Setup a Hamiltonian for the Simple Orthorhombic lattice
 *
 * This function sets up the Hamiltonian for a simple orthorhombic lattice.
 * It handles the following:
 * - Computes super-cell shape and sites
 * - Stores Hamiltonian parameters
 * - Sets local spin flags
 * - Computes transfer and interaction terms
 * - Handles different model types (spin, Hubbard, Kondo)
 * - Sets up nearest neighbor, next-nearest neighbor and third-nearest neighbor interactions
 *
 * @param[in,out] StdI Pointer to the StdIntList structure containing model parameters
 * @author Mitsuaki Kawamura (The University of Tokyo)
 */
void StdFace_Orthorhombic(
  struct StdIntList *StdI
)
{
  int isite, jsite, ntransMax, nintrMax;
  int iL, iW, iH, kCell;
  FILE *fp;
  double complex Cphase;
  double dR[3];

  /**
   * @brief Step 1: Compute super-cell shape and sites
   *
   * - Opens lattice visualization file
   * - Sets number of sites per unit cell
   * - Prints lattice parameters and dimensions
   * - Initializes site coordinates
   */
  fp = fopen("lattice.xsf", "w");
  /**/
  StdI->NsiteUC = 1;
  /**/
  fprintf(stdout, "  @ Lattice Size & Shape\n\n");
  
  StdFace_PrintVal_d("a", &StdI->a, 1.0);
  StdFace_PrintVal_d("Wlength", &StdI->length[0], StdI->a);
  StdFace_PrintVal_d("Llength", &StdI->length[1], StdI->a);
  StdFace_PrintVal_d("Hlength", &StdI->length[2], StdI->a);
  StdFace_PrintVal_d("Wx", &StdI->direct[0][0], StdI->length[0]);
  StdFace_PrintVal_d("Wy", &StdI->direct[0][1], 0.0);
  StdFace_PrintVal_d("Wz", &StdI->direct[0][2], 0.0);
  StdFace_PrintVal_d("Lx", &StdI->direct[1][0], 0.0);
  StdFace_PrintVal_d("Ly", &StdI->direct[1][1], StdI->length[1]);
  StdFace_PrintVal_d("Lz", &StdI->direct[1][2], 0.0);
  StdFace_PrintVal_d("Hx", &StdI->direct[2][0], 0.0);
  StdFace_PrintVal_d("Hy", &StdI->direct[2][1], 0.0);
  StdFace_PrintVal_d("Hz", &StdI->direct[2][2], StdI->length[2]);

  StdFace_PrintVal_d("phase0", &StdI->phase[0], 0.0);
  StdFace_PrintVal_d("phase1", &StdI->phase[1], 0.0);
  StdFace_PrintVal_d("phase2", &StdI->phase[2], 0.0);
  /**/
  StdFace_InitSite(StdI, fp, 3);
  StdI->tau[0][0] = 0.0; StdI->tau[0][1] = 0.0; ; StdI->tau[0][2] = 0.0;

  /**
   * @brief Step 2: Check and store Hamiltonian parameters
   *
   * Handles different parameter sets depending on model type:
   * - Spin model: J couplings, magnetic field, anisotropy
   * - Hubbard model: hopping terms, Coulomb interactions
   * - Kondo model: combination of both
   */
  fprintf(stdout, "\n  @ Hamiltonian \n\n");
  StdFace_NotUsed_d("K", StdI->K);
  StdFace_PrintVal_d("h", &StdI->h, 0.0);
  StdFace_PrintVal_d("Gamma", &StdI->Gamma, 0.0);
  StdFace_PrintVal_d("Gamma_y", &StdI->Gamma_y, 0.0);
  /**/
  if (strcmp(StdI->model, "spin") == 0 ) {
    StdFace_PrintVal_i("2S", &StdI->S2, 1);
    StdFace_PrintVal_d("D", &StdI->D[2][2], 0.0);
    StdFace_InputSpinNN(StdI->J, StdI->JAll, StdI->J0, StdI->J0All, "J0");
    StdFace_InputSpinNN(StdI->J, StdI->JAll, StdI->J1, StdI->J1All, "J1");
    StdFace_InputSpinNN(StdI->J, StdI->JAll, StdI->J2, StdI->J2All, "J2");
    StdFace_InputSpinNN(StdI->Jp, StdI->JpAll, StdI->J0p, StdI->J0pAll, "J0'");
    StdFace_InputSpinNN(StdI->Jp, StdI->JpAll, StdI->J1p, StdI->J1pAll, "J1'");
    StdFace_InputSpinNN(StdI->Jp, StdI->JpAll, StdI->J2p, StdI->J2pAll, "J2'");
    StdFace_InputSpin(StdI->Jpp, StdI->JppAll, "J''");
    /**/
    StdFace_NotUsed_d("mu", StdI->mu);
    StdFace_NotUsed_d("U", StdI->U);
    StdFace_NotUsed_c("t", StdI->t);
    StdFace_NotUsed_c("t0", StdI->t0);
    StdFace_NotUsed_c("t1", StdI->t1);
    StdFace_NotUsed_c("t2", StdI->t2);
    StdFace_NotUsed_c("t'", StdI->tp);
    StdFace_NotUsed_c("t0'", StdI->t0p);
    StdFace_NotUsed_c("t1'", StdI->t1p);
    StdFace_NotUsed_c("t2'", StdI->t2p);
    StdFace_NotUsed_c("t''", StdI->tpp);
    StdFace_NotUsed_d("V", StdI->V);
    StdFace_NotUsed_d("V0", StdI->V0);
    StdFace_NotUsed_d("V1", StdI->V1);
    StdFace_NotUsed_d("V'", StdI->Vp);
  }/*if (strcmp(StdI->model, "spin") == 0 )*/
  else {
    StdFace_PrintVal_d("mu", &StdI->mu, 0.0);
    StdFace_PrintVal_d("U", &StdI->U, 0.0);
    StdFace_InputHopp(StdI->t, &StdI->t0, "t0");
    StdFace_InputHopp(StdI->t, &StdI->t1, "t1");
    StdFace_InputHopp(StdI->t, &StdI->t2, "t2");
    StdFace_InputHopp(StdI->tp, &StdI->t0p, "t0'");
    StdFace_InputHopp(StdI->tp, &StdI->t1p, "t1'");
    StdFace_InputHopp(StdI->tp, &StdI->t2p, "t2'");
    StdFace_PrintVal_c("tpp", &StdI->tpp, 0.0);
    StdFace_InputCoulombV(StdI->V, &StdI->V0, "V0");
    StdFace_InputCoulombV(StdI->V, &StdI->V1, "V1");
    StdFace_InputCoulombV(StdI->V, &StdI->V2, "V2");
    StdFace_InputCoulombV(StdI->Vp, &StdI->V0p, "V0'");
    StdFace_InputCoulombV(StdI->Vp, &StdI->V1p, "V1'");
    StdFace_InputCoulombV(StdI->Vp, &StdI->V2p, "V2'");   
    StdFace_PrintVal_d("Vpp", &StdI->Vpp, 0.0);
    /**/
    StdFace_NotUsed_J("J0", StdI->J0All, StdI->J0);
    StdFace_NotUsed_J("J1", StdI->J1All, StdI->J1);
    StdFace_NotUsed_J("J2", StdI->J2All, StdI->J2);
    StdFace_NotUsed_J("J0'", StdI->J0pAll, StdI->J0p);
    StdFace_NotUsed_J("J1'", StdI->J1pAll, StdI->J1p);
    StdFace_NotUsed_J("J2'", StdI->J2pAll, StdI->J2p);
    StdFace_NotUsed_J("J0''", StdI->J0ppAll, StdI->J0pp);
    StdFace_NotUsed_J("J1''", StdI->J1ppAll, StdI->J1pp);
    StdFace_NotUsed_J("J2''", StdI->J2ppAll, StdI->J2pp);
    StdFace_NotUsed_d("D", StdI->D[2][2]);

    if (strcmp(StdI->model, "hubbard") == 0 ) {
      StdFace_NotUsed_i("2S", StdI->S2);
      StdFace_NotUsed_J("J", StdI->JAll, StdI->J);
    }/*if (strcmp(StdI->model, "hubbard") == 0 )*/
    else {
      StdFace_PrintVal_i("2S", &StdI->S2, 1);
      StdFace_InputSpin(StdI->J, StdI->JAll, "J");
    }/*if (model != "hubbard")*/
 
  }/*if (model != "spin")*/
  fprintf(stdout, "\n  @ Numerical conditions\n\n");

  /**
   * @brief Step 3: Set local spin flags and number of sites
   *
   * - Calculates total number of sites
   * - Allocates and initializes local spin flags array
   * - Handles different models (spin, Hubbard, Kondo)
   */
  StdI->nsite = StdI->NsiteUC * StdI->NCell;
  if (strcmp(StdI->model, "kondo") == 0 ) StdI->nsite *= 2;
  StdI->locspinflag = (int *)malloc(sizeof(int) * StdI->nsite);
  /**/
  if(strcmp(StdI->model, "spin") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = StdI->S2;
  else if(strcmp(StdI->model, "hubbard") == 0 )
    for (isite = 0; isite < StdI->nsite; isite++) StdI->locspinflag[isite] = 0;
  else 
    for (iL = 0; iL < StdI->nsite / 2; iL++) {
      StdI->locspinflag[iL] = StdI->S2;
      StdI->locspinflag[iL + StdI->nsite / 2] = 0;
    }

  /**
   * @brief Step 4: Calculate maximum number of interactions
   *
   * Computes upper limits for:
   * - Number of transfer terms
   * - Number of interaction terms
   * Then allocates memory for interactions
   */
  if (strcmp(StdI->model, "spin") == 0 ) {
    ntransMax = StdI->nsite * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
    nintrMax = StdI->NCell * (StdI->NsiteUC/*D*/ + 3/*J*/ + 6/*J'*/ + 4/*J''*/)
      * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
  }
  else {
    ntransMax = StdI->NCell * 2/*spin*/ * (2 * StdI->NsiteUC/*mu+h+Gamma*/ + 6/*t*/ + 12/*t'*/ + 8/*t''*/);
    nintrMax = StdI->NCell * (StdI->NsiteUC/*U*/ + 4 * (3/*V*/ + 6/*V'*/ + 4/*V''*/));

    if (strcmp(StdI->model, "kondo") == 0) {
      ntransMax += StdI->nsite / 2 * (StdI->S2 + 1/*h*/ + 2 * StdI->S2/*Gamma*/);
      nintrMax += StdI->nsite / 2 * (3 * StdI->S2 + 1) * (3 * StdI->S2 + 1);
    }/*if (strcmp(StdI->model, "kondo") == 0)*/
  }
  /**/
  StdFace_MallocInteractions(StdI, ntransMax, nintrMax);

  /**
   * @brief Step 5: Set up all interactions
   *
   * Loops through cells setting up:
   * - Local terms
   * - Nearest neighbor interactions
   * - Next-nearest neighbor interactions  
   * - Third-nearest neighbor interactions
   */
  for (kCell = 0; kCell < StdI->NCell; kCell++){
    /**/
    iW = StdI->Cell[kCell][0];
    iL = StdI->Cell[kCell][1];
    iH = StdI->Cell[kCell][2];

    /* Local terms */
    isite = kCell;
    if (strcmp(StdI->model, "kondo") == 0 ) isite += StdI->NCell;
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, -StdI->Gamma_y, isite);
      StdFace_GeneralJ(StdI, StdI->D, StdI->S2, StdI->S2, isite, isite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_HubbardLocal(StdI, StdI->mu, -StdI->h, -StdI->Gamma, -StdI->Gamma_y, StdI->U, isite);
      if (strcmp(StdI->model, "kondo") == 0 ) {
        jsite = kCell;
        StdFace_GeneralJ(StdI, StdI->J, 1, StdI->S2, isite, jsite);
        StdFace_MagField(StdI, StdI->S2, -StdI->h, -StdI->Gamma, -StdI->Gamma_y, jsite);
      }/*if (strcmp(StdI->model, "kondo") == 0 )*/
    }

    /* Nearest neighbor along W */
    StdFace_FindSite(StdI, iW, iL, iH, 1, 0, 0, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t0, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V0, isite, jsite);
    }

    /* Nearest neighbor along L */
    StdFace_FindSite(StdI, iW, iL, iH, 0, 1, 0, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J1, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t1, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V1, isite, jsite);
    }

    /* Nearest neighbor along H */
    StdFace_FindSite(StdI, iW, iL, iH, 0, 0, 1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J2, StdI->S2, StdI->S2, isite, jsite);
    }
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t2, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V2, isite, jsite);
    }

    /* Second nearest neighbor along +L+H */
    StdFace_FindSite(StdI, iW, iL, iH, 0, 1, 1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0 ) {
      StdFace_GeneralJ(StdI, StdI->J0p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t0p, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V0p, isite, jsite);
    }

    /* Second nearest neighbor along +L-H */
    StdFace_FindSite(StdI, iW, iL, iH, 0, 1, -1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J0p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t0p, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V0p, isite, jsite);
    }

    /* Second nearest neighbor along +H+W */
    StdFace_FindSite(StdI, iW, iL, iH, 1, 0, 1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J1p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t1p, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V1p, isite, jsite);
    }

    /* Second nearest neighbor along +H-W */
    StdFace_FindSite(StdI, iW, iL, iH, -1, 0, 1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J1p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t1p, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V1p, isite, jsite);
    }

    /* Second nearest neighbor along +W+L */
    StdFace_FindSite(StdI, iW, iL, iH, 1, 1, 0, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J2p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t2p, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V2p, isite, jsite);
    }

    /* Second nearest neighbor along +W-L */
    StdFace_FindSite(StdI, iW, iL, iH, 1, -1, 0, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->J2p, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->t2p, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->V2p, isite, jsite);
    }

    /* Third nearest neighbor along +W+L+H */
    StdFace_FindSite(StdI, iW, iL, iH, 1, 1, 1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->Jpp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->tpp, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->Vpp, isite, jsite);
    }

    /* Third nearest neighbor along -W+L+H */
    StdFace_FindSite(StdI, iW, iL, iH, -1, 1, 1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->Jpp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->tpp, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->Vpp, isite, jsite);
    }

    /* Third nearest neighbor along +W-L+H */
    StdFace_FindSite(StdI, iW, iL, iH, 1, -1, 1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->Jpp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->tpp, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->Vpp, isite, jsite);
    }

    /* Third nearest neighbor along +W+L-H */
    StdFace_FindSite(StdI, iW, iL, iH, 1, 1, -1, 0, 0, &isite, &jsite, &Cphase, dR);
    /**/
    if (strcmp(StdI->model, "spin") == 0) {
      StdFace_GeneralJ(StdI, StdI->Jpp, StdI->S2, StdI->S2, isite, jsite);
    }/*if (strcmp(StdI->model, "spin") == 0 )*/
    else {
      StdFace_Hopping(StdI, Cphase * StdI->tpp, isite, jsite, dR);
      StdFace_Coulomb(StdI, StdI->Vpp, isite, jsite);
    }
  }/*for (kCell = 0; kCell < StdI->NCell; kCell++)*/

  fclose(fp);
  StdFace_PrintXSF(StdI);
  StdFace_PrintGeometry(StdI);
}
