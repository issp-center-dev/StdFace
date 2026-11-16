#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "StdFace_vals.h"
#include "StdFace_ModelUtil.h"

static void release(struct StdIntList *s, int transfers, int interactions) {
  int i;
  for (i = 0; i < transfers; i++) free(s->transindx[i]);
  free(s->transindx);
  free(s->trans);
  for (i = 0; i < interactions; i++) {
    free(s->intrindx[i]);
    free(s->CintraIndx[i]);
    free(s->CinterIndx[i]);
    free(s->HundIndx[i]);
    free(s->ExIndx[i]);
    free(s->PLIndx[i]);
    free(s->PHIndx[i]);
  }
  free(s->intrindx);
  free(s->intr);
  free(s->CintraIndx);
  free(s->Cintra);
  free(s->CinterIndx);
  free(s->Cinter);
  free(s->HundIndx);
  free(s->Hund);
  free(s->ExIndx);
  free(s->Ex);
  free(s->PLIndx);
  free(s->PairLift);
  free(s->PHIndx);
  free(s->PairHopp);
}

int main(int argc, char **argv) {
  struct StdIntList s;
  int model, mask, site, i, failed = 0;
  if (argc > 1) {
    memset(&s, 0, sizeof(s));
    strcpy(s.model, "hubbard");
    s.nsite = 2;
    if (strcmp(argv[1], "negative-transfers") == 0) {
      StdFace_MallocInteractions(&s, -1, 2);
    } else if (strcmp(argv[1], "negative-interactions") == 0) {
      StdFace_MallocInteractions(&s, 8, -1);
    } else if (strcmp(argv[1], "negative-sites") == 0) {
      s.nsite = -1;
      StdFace_MallocInteractions(&s, 8, 2);
    } else if (strcmp(argv[1], "overflow") == 0) {
      StdFace_MallocInteractions(&s, INT_MAX - 1, 2);
    } else {
      return 2;
    }
    return 0; /* Invalid inputs must exit before returning here. */
  }
  for (model = 0; model < 2; model++) {
    /* Start with all fields nonzero to exercise the formerly undersized buffer. */
    for (mask = 7; mask >= 0; mask--) {
      double h = (mask & 1) ? 0.2 : 0;
      double gx = (mask & 2) ? 0.4 : 0;
      double gy = (mask & 4) ? 0.6 : 0;
      double complex actual[3][2][2] = {0};
      memset(&s, 0, sizeof(s));
      strcpy(s.model, model ? "kondo" : "hubbard");
      s.nsite = model ? 6 : 3;
      StdFace_MallocInteractions(&s, 4 * 3, 3);
      for (site = 0; site < 3; site++) {
        StdFace_HubbardLocal(&s, 0.31, h, gx, gy, 4, site);
      }
      for (i = 0; i < s.ntrans; i++) {
        int *p = s.transindx[i];
        if (p[0] < 0 || p[0] >= 3 || p[0] != p[2] ||
            p[1] < 0 || p[1] > 1 || p[3] < 0 || p[3] > 1) {
          failed = 1;
          break;
        }
        actual[p[0]][p[1]][p[3]] += s.trans[i];
      }
      /* Independent spin-1/2 matrix, summed over separate real/imaginary entries. */
      for (site = 0; site < 3; site++) {
        failed |= cabs(actual[site][0][0] - (0.31 - 0.5 * h)) > 1e-14;
        failed |= cabs(actual[site][1][1] - (0.31 + 0.5 * h)) > 1e-14;
        failed |= cabs(actual[site][1][0] - (-0.5 * gx - 0.5 * I * gy)) > 1e-14;
        failed |= cabs(actual[site][0][1] - (-0.5 * gx + 0.5 * I * gy)) > 1e-14;
      }
      failed |= s.NCintra != 3;
      release(&s, 6 * 3, 3);
    }
  }
  if (failed) {
    fprintf(stderr, "Local transfer coefficients differ from the independent spin matrices.\n");
  }
  return failed;
}
