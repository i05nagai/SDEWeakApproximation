/*
 * sde_wa.h
 */
/*
 * $Source$
 * $Revision$
 * $Author$
 * $Date$
 */

#ifndef _SDE_WA_H_
#define _SDE_WA_H_

#include "sde_wa_errno.h"

enum SDE_type {ITO=0, STR=1};/* Ito=>0, Stratonovich=>1 */

typedef struct {
  enum SDE_type sde_type;
  int (**V)(const double y[], double dy[], void *params);
  int (**drift_corrector)(const double y[], double *dVdy[], void *params);
  int (**exp_sV)(double s, const double y[], double exp_sVy[], void *params);
  int dim_y;
  int dim_BM;
  void *params;
} SDE_WA_SYSTEM;

enum ALG {E_M=0, N_V=1, N_N=2};
enum Bernoulli_rv {T=0, H=1};

typedef struct {
  enum Bernoulli_rv rv_nv_b;
  double *rv_nv_n;
} RV_NV;



typedef struct sde_wa_sltn{
  enum ALG alg;
  int mth_is;
  SDE_WA_SYSTEM *sde;
  int (*one_step)(struct sde_wa_sltn *sl, double s);
  double *initv;
  double *destv;
  double *drift_step_interv; /* intermidiate values in drift correction etc. */
  double **drift_corrector_matrix; /* for Ito<->Stratonovich */
  union{
    double *em; 
    RV_NV *nv; 
    double *nn;
  } sample_pt;
  double *rk_step_interv; /* RK intermidiate values for NV or NN*/
  double *nn_sample_pt_interv;
} SDE_WA_SLTN;

SDE_WA_SYSTEM *alloc_SDE_WA_SYSTEM(int N, int d, void *params);
void free_SDE_WA_SYSTEM(SDE_WA_SYSTEM *sys);
SDE_WA_SLTN *alloc_SDE_WA_SLTN(enum ALG alg, int mth_is, SDE_WA_SYSTEM *sde);
void free_SDE_WA_SLTN(SDE_WA_SLTN *sltn);

void Ito_to_Strt_drift(SDE_WA_SLTN *X, const double *init, double *dest);

int next_SDE_WA(SDE_WA_SLTN *X, double s, double y[], double dy[], void *rv);

#endif
