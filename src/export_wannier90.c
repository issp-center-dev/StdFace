#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "StdFace_vals.h"
#include "StdFace_ModelUtil.h"

#if defined(_HWAVE)

// #undef NDEBUG

#define fatal(fmt,...)                                    \
  do {                                                    \
    fprintf(stderr, "ERROR: %s: " fmt " in %s:%d\n",      \
            __func__, ##__VA_ARGS__, __FILE__, __LINE__); \
    StdFace_exit(-1);                                     \
  } while(0)


/* control flags */
int is_export_all = 1;  /* default for exportall parameter */

/* parameters */
static const double eps = 1.0e-8;
  
/* type declarations */
typedef struct {
  int r[3];   /* relative coordinate */
  int a, b;   /* orbit index */
  int s, t;   /* extension: spin dependence */
  double complex v;   /* value */
} IntrItem;


/**
   @brief Print geometry data in Wannier90 geom format
*/
void WriteGeometry(struct StdIntList *StdI, char *fname)
{
  FILE *fp = fopen(fname, "w");
  if (!fp) fatal("cannot open file for output: %s", fname);

  /* print primitive vectors */
  for (int ii = 0; ii < 3; ++ii) {
    fprintf(fp, "%16.12f %16.12f %16.12f\n",
            StdI->direct[ii][0], StdI->direct[ii][1], StdI->direct[ii][2]);
  }

  /* print number of orbits */
  fprintf(fp, "%d\n", StdI->NsiteUC);

  /* print center of orbits */
  for (int k = 0; k < StdI->NsiteUC; ++k) {
    fprintf(fp, "%25.15e %25.15e %25.15e\n",
            StdI->tau[k][0], StdI->tau[k][1], StdI->tau[k][2]);
  }

  fflush(fp);
  fclose(fp);

  printf("%24s is written.\n", fname);

  return;
}

/**
   @brief Print interaction coefficient data in Wannier90 (-like) format
*/
void WriteWannier90(int nintr_table, IntrItem *intr_table,
                    int nsiteuc, int nspin,
                    char *fname, char *tagname)

{
  /* find range */
  int rmin[3], rmax[3];
  for (int i = 0; i < 3; ++i) {
    rmin[i] = intr_table[0].r[i];
    rmax[i] = intr_table[0].r[i];
  }

  for (int k = 1; k < nintr_table; ++k) {
    for (int i = 0; i < 3; ++i) {
      int r = intr_table[k].r[i];
      if (r < rmin[i]) rmin[i] = r;
      if (r > rmax[i]) rmax[i] = r;
    }
  }

  int rr[3];
  for (int i = 0; i < 3; ++i) {
    rr[i] = (abs(rmin[i]) > abs(rmax[i])) ? abs(rmin[i]) : abs(rmax[i]);
  }
  
#ifndef NDEBUG
  for (int i = 0; i < 3; ++i) {
    printf("DEBUG: %s: range_%d = [%d, %d], rr=%d\n", __func__, i, rmin[i], rmax[i], rr[i]);
  }
#endif

  int nvol = (rr[0]*2+1) * (rr[1]*2+1) * (rr[2]*2+1);
  int matrix_size = nvol * nsiteuc * nsiteuc * nspin * nspin;

  double complex *matrix = (double complex *)malloc(sizeof(double complex) * matrix_size);
  if (!matrix) fatal("cannot allocate matrix");

  for (int i = 0; i < matrix_size; ++i) matrix[i] = 0.0;


#define _index(rx,ry,rz,a,b,s,t)                \
  ((t) + nspin * (                              \
    (s) + nspin * (                             \
      (b) + nsiteuc * (                         \
        (a) + nsiteuc * (                       \
          ((rz)+rr[2]) + (rr[2]*2+1) * (        \
            ((ry)+rr[1]) + (rr[1]*2+1) * (      \
              ((rx)+rr[0])                      \
              )))))))

  for (int k = 0; k < nintr_table; ++k) {

    matrix[_index( intr_table[k].r[0],
                   intr_table[k].r[1],
                   intr_table[k].r[2],
                   intr_table[k].a,
                   intr_table[k].b,
                   intr_table[k].s,
                   intr_table[k].t)] = intr_table[k].v;

    int rindex = _index(-intr_table[k].r[0],
                        -intr_table[k].r[1],
                        -intr_table[k].r[2],
                        intr_table[k].b,
                        intr_table[k].a,
                        intr_table[k].t,
                        intr_table[k].s);

    if (cabs(matrix[rindex]) < eps) {
      matrix[rindex] = conj(intr_table[k].v);
    }
  }

#ifndef NDEBUG
  /* dump */
  for (int r = 0; r < nvol; ++r) {
    int rz = (r) % (rr[2]*2+1) - rr[2];
    int ry = (r / (rr[2]*2+1)) % (rr[1]*2+1) - rr[1];
    int rx = (r / ((rr[2]*2+1)*(rr[1]*2+1))) % (rr[0]*2+1) - rr[0];

    for (int a = 0; a < nsiteuc; ++a) {
      for (int b = 0; b < nsiteuc; ++b) {
        for (int s = 0; s < nspin; ++s) {
          for (int t = 0; t < nspin; ++t) {
            printf("DEBUG: J_abst(r) : r=[%2d,%2d,%2d], a,b=%2d,%2d, s,t=%2d,%2d : %f %f\n",
                   rx, ry, rz, a, b, s, t,
                   creal(matrix[_index(rx,ry,rz,a,b,s,t)]),
                   cimag(matrix[_index(rx,ry,rz,a,b,s,t)]));
          }
        }
      }
    }
  }
#endif

  FILE *fp = fopen(fname, "w");
  if (!fp) fatal("cannot open file: %s", fname);

  /* write header */
  fprintf(fp, "%s in wannier90-like format for uhfk\n", tagname);
  fprintf(fp, "%d\n%d\n", nsiteuc, nvol);
  for (int i = 0; i < nvol; ++i) {
    if (i > 0 && i % 15 == 0) fprintf(fp, "\n");
    fprintf(fp, " %d", 1);
  }
  fprintf(fp, "\n");

  /* write matrix body */
  for (int r = 0; r < nvol; ++r) {
    int rz = (r) % (rr[2]*2+1) - rr[2];
    int ry = (r / (rr[2]*2+1)) % (rr[1]*2+1) - rr[1];
    int rx = (r / ((rr[2]*2+1)*(rr[1]*2+1))) % (rr[0]*2+1) - rr[0];

    for (int a = 0; a < nsiteuc; ++a) {
      for (int b = 0; b < nsiteuc; ++b) {

        if (nspin > 1) {
          /* extended format */
          for (int s = 0; s < nspin; ++s) {
            for (int t = 0; t < nspin; ++t) {
              if (is_export_all || cabs(matrix[_index(rx,ry,rz,a,b,s,t)]) > eps) {
                /* note: orbit a,b starts 1 */
                fprintf(fp, "%4d %4d %4d %4d %4d %4d %4d %16.12f %16.12f\n",
                        rx, ry, rz, (a+1), (b+1), s, t, 
                        creal(matrix[_index(rx,ry,rz,a,b,s,t)]),
                        cimag(matrix[_index(rx,ry,rz,a,b,s,t)]));
              }
            }
          }
        } else {
          int s = 0;
          int t = 0;
          if (is_export_all || cabs(matrix[_index(rx,ry,rz,a,b,s,t)]) > eps) {
            /* note: orbit a,b starts 1 */
            fprintf(fp, "%4d %4d %4d %4d %4d %16.12f %16.12f\n",
                    rx, ry, rz, (a+1), (b+1),
                    creal(matrix[_index(rx,ry,rz,a,b,s,t)]),
                    cimag(matrix[_index(rx,ry,rz,a,b,s,t)]));
          }
        }

      }
    }
  }

  fflush(fp);
  fclose(fp);

  printf("%24s is written.\n", fname);

#undef _index

  /* tidy up */
  free(matrix);

  return;
}

/**
   @bridf unfold site from [0, N] to [-N/2, N/2]
*/
void unfold_site(struct StdIntList *StdI, int v_in[3], int v_out[3])
{
  int v[3];
  for (int i = 0; i < 3; ++i) {
    v[i] = 0;
    for (int j = 0; j < 3; ++j) {
      v[i] += StdI->rbox[i][j] * v_in[j];
    }
  }

  for (int i = 0; i < 3; ++i) {
    double vv = 1.0 * v[i] / StdI->NCell;
    if (vv > 0.5) {
      v[i] -= StdI->NCell;
    } else if (vv <= -0.5) {
      v[i] += StdI->NCell;
    }
  }

  int w[3];
  for (int i = 0; i < 3; ++i) {
    w[i] = 0;
    for (int j = 0; j < 3; ++j) {
      w[i] += v[j] * StdI->box[j][i];
    }
    w[i] /= StdI->NCell;
  }

  for (int i = 0; i < 3; ++i) {
    v_out[i] = w[i];
  }

  return;
}

/**
   @brief generate key as a pair for interaction table
   ordered = 0 -> as-is
   ordered = 1 -> i < j
*/
void generate_key(int keylen, int *index_out, int *index, int ordered)
{
  if (keylen == 1) {  /* i */
    index_out[0] = index[0];

  } else if (keylen == 2) {  /* i, j */
    int i = index[0];
    int j = index[1];
    if (ordered == 1 && i > j) {
      index_out[0] = j;
      index_out[1] = i;
    } else {
      index_out[0] = i;
      index_out[1] = j;
    }

  } else if (keylen == 4) {  /* (i,s), (j,t) */
    int i = index[0];
    int s = index[1];
    int j = index[2];
    int t = index[3];
    if (ordered == 1 && i > j) {
      index_out[0] = j;
      index_out[1] = t;
      index_out[2] = i;
      index_out[3] = s;
    } else {
      index_out[0] = i;
      index_out[1] = s;
      index_out[2] = j;
      index_out[3] = t;
    }

  } else {
    fatal("unsupported keylen: %d", keylen);
  }
  return;
}

/**
   @brief copy key data
*/
void store_key(int keylen, int *dst_key, int *src_key)
{
  for (int k = 0; k < keylen; ++k) {
    dst_key[k] = src_key[k];
  }
  return;
}

/**
   @brief test whether two keys are same
*/
int is_equal_key(int keylen, int *key_a, int *key_b)
{
  for (int k = 0; k < keylen; ++k) {
    if (key_a[k] != key_b[k]) return 0;
  }
  return 1;
}

/**
   @brief make string representation of a key
*/
char *to_string_key(int keylen, int *key)
{
  static char buf[1024];

  if (keylen == 1) {
    snprintf(buf, 1024, "(%d)", key[0]);
  } else if (keylen == 2) {
    snprintf(buf, 1024, "(%d, %d)", key[0], key[1]);
  } else if (keylen == 4) {
    snprintf(buf, 1024, "(%d, %d, %d, %d)", key[0], key[1], key[2], key[3]);
  } else {
    fatal("unsupported keylen: %d", keylen);
  }

  return buf;
}

/**
   @brief make table compact by accumulating entries of the same key
*/
void accumulate_list(int keylen,
                     int ntbl, int **tbl_index, double complex *tbl_value,
                     int *ntbl_out, int **tbl_index_out, double complex *tbl_value_out,
                     int ordered)
{
  *ntbl_out = 0;

  for (int k = 0; k < ntbl; ++k) {

    int idx[keylen];
    generate_key(keylen, idx, tbl_index[k], ordered);
    /* i.e. idx = tbl_index[k]; */
    double complex val = tbl_value[k];

#ifndef NDEBUG
    printf("DEBUG: %s: k=%d, idx=%s, v=(%f,%f)\n", __func__, k, to_string_key(keylen, idx), creal(val), cimag(val));
#endif
    
    int is_found = 0;
    int jj = 0;
    for (int j = 0; j < *ntbl_out; ++j) {
      if (is_equal_key(keylen, tbl_index_out[j], idx)) {
        is_found = 1;
        jj = j;
        break;
      }
    }

    if (is_found) {
      tbl_value_out[jj] += val;
#ifndef NDEBUG
      printf("DEBUG: %s: key found.\n", __func__);
#endif
    } else {
      store_key(keylen, tbl_index_out[*ntbl_out], idx);
      tbl_value_out[*ntbl_out] = val;
      *ntbl_out += 1;
#ifndef NDEBUG
      printf("DEBUG: %s: new entry.\n", __func__);
#endif
    }
  }

  /* eliminate zero */
  for (int k = *ntbl_out - 1; k >= 0; --k) {

    if (cabs(tbl_value_out[k]) < eps) {
#ifndef NDEBUG
      printf("DEBUG: %s: eliminate entry %d: index=%s, value=(%f,%f).\n",
             __func__,
             k,
             to_string_key(keylen, tbl_index_out[k]),
             creal(tbl_value_out[k]), cimag(tbl_value_out[k]));
#endif
      for (int j = k+1; j < *ntbl_out; ++j) {
        store_key(keylen, tbl_index_out[j-1], tbl_index_out[j]);
        tbl_value_out[j-1] = tbl_value_out[j];
      }
      *ntbl_out -= 1;
    }
  }

#ifndef NDEBUG
  /* dump */
  printf("DEBUG: %s: ntbl_out = %d\n", __func__, *ntbl_out);
  for (int k = 0; k < *ntbl_out; ++k) {
    printf("DEBUG: %s: %d: index=%s, value=(%f,%f)\n",
           __func__,
           k,
           to_string_key(keylen, tbl_index_out[k]),
           creal(tbl_value_out[k]), cimag(tbl_value_out[k]));
  }
#endif

  return;
}

/**
   @brief Print interaction coefficient of StdFace internal data structure (complex array)
*/
void ExportInter(struct StdIntList *StdI,
                 int ntbl, int **tbl_index, double complex *tbl_value,
                 char *fname, char *tagname)
{
  if (ntbl == 0) {
    printf("%24s is skipped.\n", fname);
    return;
  }

  /* accumulate entries of the same index pair */
  int *intr_index_buf = (int *)malloc(sizeof(int) * ntbl * 2);
  if (!intr_index_buf) fatal("cannot allocate intr_index_buf");

  int **intr_index = (int **)malloc(sizeof(int *) * ntbl);
  if (!intr_index) fatal("cannot allocate intr_index");

  for (int k = 0; k < ntbl; ++k) {
    intr_index[k] = intr_index_buf + 2 * k;
  }

  double complex *intr_value = (double complex *)malloc(sizeof(double complex) * ntbl);
  if (!intr_value) fatal("cannot allocate intr_value");

  int nintr = 0;

  accumulate_list(2 /* keylen */,
                  ntbl, tbl_index, tbl_value,
                  &nintr, intr_index, intr_value,
                  1); /* i < j pair */

  if (nintr > 0) {

    /* elements in relative coordinate */
    IntrItem *intr_table = (IntrItem *)malloc(sizeof(IntrItem) * nintr);
    if (!intr_table) fatal("cannot allocate intr");

    int nintr_table = 0;

    /* relative site index */
    for (int k = 0; k < nintr; ++k) {
      int idx_i = intr_index[k][0];
      int idx_j = intr_index[k][1];

      int icell = idx_i / StdI->NsiteUC;
      int isite = idx_i % StdI->NsiteUC;
      int jcell = idx_j / StdI->NsiteUC;
      int jsite = idx_j % StdI->NsiteUC;

      int rr[3];
      for (int i = 0; i < 3; ++i) {
        rr[i] = StdI->Cell[jcell][i] - StdI->Cell[icell][i];
      }
      unfold_site(StdI, rr, rr);

      /* check consistency */
      int is_found = 0;
      int jj = 0;
      for (int j = 0; j < nintr_table; ++j) {
        if (   intr_table[j].r[0] == rr[0]
               && intr_table[j].r[1] == rr[1]
               && intr_table[j].r[2] == rr[2]
               && intr_table[j].a == isite
               && intr_table[j].b == jsite
          ){
          is_found = 1;
          jj = j;

          if (cabs(intr_table[j].v - intr_value[k]) > eps) {
            printf("WARNING: not uniform. expected=(%f,%f), found=(%f,%f) for index %d,%d\n",
                   creal(intr_table[j].v), cimag(intr_table[j].v),
                   creal(intr_value[k]), cimag(intr_value[k]),
                   idx_i, idx_j);
          }
          break;
        }
      }

      /* store */
      if (!is_found) {
        for (int i = 0; i < 3; ++i) {
          intr_table[nintr_table].r[i] = rr[i];
        }
        intr_table[nintr_table].a = isite;
        intr_table[nintr_table].b = jsite;
        intr_table[nintr_table].s = 0;
        intr_table[nintr_table].t = 0;
        intr_table[nintr_table].v = intr_value[k];
        ++nintr_table;
      }
    }

#ifndef NDEBUG
    /* dump */
    printf("nintr_table=%d\n", nintr_table);
    for (int k = 0; k < nintr_table; ++k) {
      printf("DEBUG: intr_table[%d]: r=[%2d,%2d,%2d], a,b=%2d,%2d, s,t=%2d,%2d, val=%f,%f\n",
             k,
             intr_table[k].r[0], intr_table[k].r[1], intr_table[k].r[2],
             intr_table[k].a, intr_table[k].b, intr_table[k].s, intr_table[k].t,
             creal(intr_table[k].v), cimag(intr_table[k].v));
    }
#endif

    /* write to file */
    WriteWannier90(nintr_table, intr_table, StdI->NsiteUC, 1, fname, tagname);
  
    /* tidy up */
    free(intr_table);

  } else {
    printf("%24s is skipped.\n", fname);
  } /* if (nintr > 0) */
  
  free(intr_index_buf);
  free(intr_index);
  free(intr_value);

  return;
}

/**
   @brief Print interaction coefficient of StdFace internal data structure (real array)
*/
void ExportInterReal(struct StdIntList *StdI,
                     int ntbl, int **tbl_index, double *tbl_value,
                     char *fname, char *tagname)
{
  double complex *buf = (double complex *)malloc(sizeof(double complex) * ntbl);
  if (!buf) fatal("cannot allocate buffer");

  for (int k = 0; k < ntbl; ++k) {
    buf[k] = tbl_value[k];
  }

  ExportInter(StdI, ntbl, tbl_index, buf, fname, tagname);

  /* tidy up */
  free(buf);

  return;
}

/**
   @brief Print coefficient of transfer term
*/
void ExportTransfer(struct StdIntList *StdI,
                    int ntbl, int **tbl_index, double complex *tbl_value,
                    char *fname, char *tagname, int spin_dep)
{
  if (ntbl == 0) {
    printf("%24s is skipped.\n", fname);
    return;
  }

  /* accumulate entries of the same index pair */
  int *intr_index_buf = (int *)malloc(sizeof(int) * ntbl * 4);
  if (!intr_index_buf) fatal("cannot allocate intr_index_buf");

  int **intr_index = (int **)malloc(sizeof(int *) * ntbl);
  if (!intr_index) fatal("cannot allocate intr_index");

  for (int k = 0; k < ntbl; ++k) {
    intr_index[k] = intr_index_buf + 4 * k;
  }

  double complex *intr_value = (double complex *)malloc(sizeof(double complex) * ntbl);
  if (!intr_value) fatal("cannot allocate intr_value");

  int nintr = 0;

  accumulate_list(4 /* keylen */,
                  ntbl, tbl_index, tbl_value,
                  &nintr, intr_index, intr_value, 0);  /* not accumulate to i<j pair */

  if (nintr > 0) {

    /* elements in relative coordinate */
    IntrItem *intr_table = (IntrItem *)malloc(sizeof(IntrItem) * nintr);
    if (!intr_table) fatal("cannot allocate intr");

    int nintr_table = 0;

    /* relative site index */
    for (int k = 0; k < nintr; ++k) {
      int idx_i = intr_index[k][0];
      int ispin = intr_index[k][1];
      int idx_j = intr_index[k][2];
      int jspin = intr_index[k][3];

      int icell = idx_i / StdI->NsiteUC;
      int isite = idx_i % StdI->NsiteUC;
      int jcell = idx_j / StdI->NsiteUC;
      int jsite = idx_j % StdI->NsiteUC;

      int rr[3];
      for (int i = 0; i < 3; ++i) {
        rr[i] = StdI->Cell[jcell][i] - StdI->Cell[icell][i];
      }
      unfold_site(StdI, rr, rr);

      intr_value[k] *= -1;  /* by convention */

      /* check consistency */
      int is_found = 0;
      int jj = 0;
      for (int j = 0; j < nintr_table; ++j) {
        if (   intr_table[j].r[0] == rr[0]
               && intr_table[j].r[1] == rr[1]
               && intr_table[j].r[2] == rr[2]
               && intr_table[j].a == isite
               && intr_table[j].b == jsite
               && intr_table[j].s == ispin
               && intr_table[j].t == jspin
          ){
          is_found = 1;
          jj = j;

          if (cabs(intr_table[j].v - intr_value[k]) > eps) {
            printf("WARNING: not uniform. expected=(%f,%f), found=(%f,%f) for index %d,%d\n",
                   creal(intr_table[j].v), cimag(intr_table[j].v),
                   creal(intr_value[k]), cimag(intr_value[k]),
                   idx_i, idx_j);
          }
          break;
        }
      }

      /* store */
      if (!is_found) {

        if (spin_dep == 0 && !(ispin == 0 && jspin == 0)) {
          continue;  /* skip */
        }

        for (int i = 0; i < 3; ++i) {
          intr_table[nintr_table].r[i] = rr[i];
        }
        intr_table[nintr_table].a = isite;
        intr_table[nintr_table].b = jsite;
        intr_table[nintr_table].s = ispin;
        intr_table[nintr_table].t = jspin;
        intr_table[nintr_table].v = intr_value[k];
        ++nintr_table;
      }
    }

#ifndef NDEBUG
    /* dump */
    printf("nintr_table=%d\n", nintr_table);
    for (int k = 0; k < nintr_table; ++k) {
      printf("DEBUG: intr_table[%d]: r=[%2d,%2d,%2d], a,b=%2d,%2d, s,t=%2d,%2d, val=%f,%f\n",
             k,
             intr_table[k].r[0], intr_table[k].r[1], intr_table[k].r[2],
             intr_table[k].a, intr_table[k].b, intr_table[k].s, intr_table[k].t,
             creal(intr_table[k].v), cimag(intr_table[k].v));
    }
#endif

    /* write to file */
    WriteWannier90(nintr_table, intr_table, StdI->NsiteUC, ((spin_dep == 1) ? 2 : 1), fname, tagname);
  
    /* tidy up */
    free(intr_table);

  } else {
    printf("%24s is skipped.\n", fname);
  } /* if (nintr > 0) */
  
  free(intr_index_buf);
  free(intr_index);
  free(intr_value);

  return;
}

/**
   @brief Print coefficient of on-site Coulomb term
*/
void ExportCoulombIntra(struct StdIntList *StdI,
                        int ntbl, int **tbl_index, double *tbl_value,
                        char *fname, char *tagname)
{
  if (ntbl == 0) {
    printf("%24s is skipped.\n", fname);
    return;
  }

  /* convert to complex array */
  double complex *tbl_value_c = (double complex *)malloc(sizeof(double complex) * ntbl);
  if (!tbl_value_c) fatal("cannot allocate tbl_value_c");

  for (int i = 0; i < ntbl; ++i) {
    tbl_value_c[i] = tbl_value[i];
  }

  /* accumulate entries of the same index */
  int *intr_index_buf = (int *)malloc(sizeof(int) * ntbl * 1);
  if (!intr_index_buf) fatal("cannot allocate intr_index_buf");

  int **intr_index = (int **)malloc(sizeof(int *) * ntbl);
  if (!intr_index) fatal("cannot allocate intr_index");

  for (int k = 0; k < ntbl; ++k) {
    intr_index[k] = intr_index_buf + 1 * k;
  }

  double complex *intr_value = (double complex *)malloc(sizeof(double complex) * ntbl);
  if (!intr_value) fatal("cannot allocate intr_value");

  int nintr = 0;

  accumulate_list(1 /* keylen */,
                  ntbl, tbl_index, tbl_value_c,
                  &nintr, intr_index, intr_value, 1);

  /* check if uniform */
  if (nintr > 0) {

    IntrItem *intr_table = (IntrItem *)malloc(sizeof(IntrItem) * nintr);
    if (!intr_table) fatal("cannot allocate intr");

    int nintr_table = 0;

    for (int k = 0; k < nintr; ++k) {
      int idx_i = intr_index[k][0];
      int icell = idx_i / StdI->NsiteUC;
      int isite = idx_i % StdI->NsiteUC;

      int is_found = 0;
      int jj = 0;
      for (int j = 0; j < nintr_table; ++j) {
        if (intr_table[j].a == isite) {
          is_found = 1;
          jj = j;

          if (cabs(intr_table[j].v - intr_value[k]) > eps) {
            printf("WARNING: not uniform. expected=(%f,%f), found=(%f,%f) for index %d\n",
                   creal(intr_table[j].v), cimag(intr_table[j].v),
                   creal(intr_value[k]), cimag(intr_value[k]),
                   idx_i);
          }
          break;
        }
      }

      if (!is_found) {
        intr_table[nintr_table].r[0] = 0;
        intr_table[nintr_table].r[1] = 0;
        intr_table[nintr_table].r[2] = 0;
        intr_table[nintr_table].a = isite;
        intr_table[nintr_table].b = isite;
        intr_table[nintr_table].s = 0;
        intr_table[nintr_table].t = 0;
        intr_table[nintr_table].v = intr_value[k];
        ++nintr_table;
      }
    }

    WriteWannier90(nintr_table, intr_table, StdI->NsiteUC, 1, fname, tagname);

    /* tidy up */
    free(intr_table);

  } else {
    printf("%24s is skipped.\n", fname);
  }

  /* tidy up */
  free(tbl_value_c);
  free(intr_index_buf);
  free(intr_index);
  free(intr_value);

  return;
}

/**
   @brief add prefix to filename if parameter fileprefix is set.
*/
char *prefix_(struct StdIntList *StdI, char *fname)
{
  static char buf[4096];

  if (strcmp(StdI->fileprefix, "****")==0) {
    strcpy(buf, fname);
  } else {
    snprintf(buf, 4096, "%s_%s", StdI->fileprefix, fname);
  }
  return buf;
}

/**
   @brief Print geometry information
*/
void ExportGeometry(struct StdIntList *StdI)
{
  return WriteGeometry(StdI, prefix_(StdI, "geom.dat"));
}

/**
   @brief Print interaction term coefficients
*/
void ExportInteraction(struct StdIntList *StdI)
{
  if (StdI->export_all != StdI->NaN_i) is_export_all = StdI->export_all;

  ExportTransfer(StdI,
                 StdI->ntrans, StdI->transindx, StdI->trans,
                 prefix_(StdI, "transfer.dat"), "Transfer",
                 0);

  ExportCoulombIntra(StdI,
                     StdI->NCintra, StdI->CintraIndx, StdI->Cintra,
                     prefix_(StdI, "coulombintra.dat"), "CoulombIntra");

  ExportInterReal(StdI,
                  StdI->NCinter, StdI->CinterIndx, StdI->Cinter,
                  prefix_(StdI, "coulombinter.dat"), "CoulombInter");
  ExportInterReal(StdI,
                  StdI->NHund, StdI->HundIndx, StdI->Hund,
                  prefix_(StdI, "hund.dat"), "Hund");
  ExportInterReal(StdI,
                  StdI->NEx, StdI->ExIndx, StdI->Ex,
                  prefix_(StdI, "exchange.dat"), "Exchange");
  ExportInterReal(StdI,
                  StdI->NPairLift, StdI->PLIndx, StdI->PairLift,
                  prefix_(StdI, "pairlift.dat"), "PairLift");
  ExportInterReal(StdI,
                  StdI->NPairHopp, StdI->PHIndx, StdI->PairHopp,
                  prefix_(StdI, "pairhopp.dat"), "PairHopp");

  return;
}

#endif /* defined(_HWAVE) */
