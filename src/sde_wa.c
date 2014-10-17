/*
 * sde_wa.c
 */
/*
 * $Source$
 * $Revision$
 * $Author$
 * $Date$
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; version 2.1 of the License.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License version 2.1 for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * version 2.1 along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, US
 */

#include <stdlib.h>

#include "sde_wa.h"
#include "sde_wa_em.h"
#include "sde_wa_nv.h"
#include "sde_wa_nn.h"

SDE_WA_SYSTEM *alloc_SDE_WA_SYSTEM(int N, int d, void *par){
  SDE_WA_SYSTEM *sde;
  
  sde=(SDE_WA_SYSTEM *)malloc(sizeof(SDE_WA_SYSTEM));
  sde->V
    =(int (**)(const double y[], double dy[], void *params))
    malloc(sizeof(int (*)(const double y[], double dy[], void *params))*(1+d));
  sde->drift_corrector
    =(int (**)(const double y[], double *dVdy[], void *params))
    malloc(sizeof(int (*)(const double y[],
			  double *dVdy[],
			  void *params))*(1+d));
  sde->exp_sV
    =(int (**)(double s, const double y[], double exp_sVy[], void *params))
    malloc(sizeof(int (*)(double s,
			  const double y[],
			  double exp_sVy[],
			  void *params))*(1+d));
  sde->dim_y=N;
  sde->dim_BM=d;
  sde->params=par;
  return sde;
}

void free_SDE_WA_SYSTEM(SDE_WA_SYSTEM *sde){
  free(sde->V);
  free(sde->drift_corrector);
  free(sde->exp_sV);
  free(sde);
}

SDE_WA_SLTN *alloc_SDE_WA_SLTN(enum ALG alg, int mth_is, SDE_WA_SYSTEM *sde){
  
  int j;
  SDE_WA_SLTN *sltn;
  sltn=(SDE_WA_SLTN *)malloc(sizeof(SDE_WA_SLTN));

  sltn->mth_is=mth_is;
  sltn->alg=alg;
  sltn->sde=sde;
  sltn->initv=(double *)malloc(sizeof(double)*sde->dim_y);
  sltn->destv=(double *)malloc(sizeof(double)*sde->dim_y);
  sltn->drift_corrector_matrix=(double **)malloc(sizeof(double *)*sde->dim_y);
  sltn->drift_step_interv=(double *)malloc(sizeof(double)*(sde->dim_y));
  for (j=0; j<sde->dim_y; j++)
    sltn->drift_corrector_matrix[j]=(double *)malloc(sizeof(double)*sde->dim_y);
  switch (alg){
  case 0:
    sltn->one_step=sde_wa_em;
    sltn->sample_pt.em=(double *)malloc(sizeof(double)*sde->dim_BM);
    break;
  case 1:
    sltn->one_step=sde_wa_nv;
    if (sltn->mth_is==5)
      sltn->rk_step_interv=(double *)malloc(sizeof(double)*(sde->dim_y)*11);
    else
      sltn->rk_step_interv=(double *)malloc(sizeof(double)*(sde->dim_y)*17);
    sltn->sample_pt.nv=(RV_NV *)malloc(sizeof(RV_NV));
    sltn->sample_pt.nv->rv_nv_n=(double *)malloc(sizeof(double)*(sde->dim_BM));
    break;
  case 2:
    sltn->one_step=sde_wa_nn;
    if (sltn->mth_is==5)
      sltn->rk_step_interv=(double *)malloc(sizeof(double)*(sde->dim_y)*11);
    else
      sltn->rk_step_interv=(double *)malloc(sizeof(double)*(sde->dim_y)*17);
    sltn->sample_pt.nn=(double *)malloc(sizeof(double)*(2*sde->dim_BM));
    sltn->nn_sample_pt_interv=(double *)malloc(sizeof(double)*(2*sde->dim_BM));
    break;
  }
  return sltn;
}

void free_SDE_WA_SLTN(SDE_WA_SLTN *sltn){
  int j;
  free(sltn->initv);
  free(sltn->drift_step_interv);
  free(sltn->destv);
  for (j=0; j<sltn->sde->dim_y;j++) free(sltn->drift_corrector_matrix[j]);
  free(sltn->drift_corrector_matrix);
  switch (sltn->alg){
  case 0:
    free(sltn->sample_pt.em);
    break;
  case 1:
    free(sltn->rk_step_interv);
    free(sltn->sample_pt.nv->rv_nv_n);
    free(sltn->sample_pt.nv);
    break;
  case 2:
    free(sltn->rk_step_interv);
    free(sltn->sample_pt.nn);
    free(sltn->nn_sample_pt_interv);
    break;
  }
  free(sltn);
}

void Ito_to_Strt_drift(SDE_WA_SLTN *X, const double *init, double *dest){
  
  int i,j,k;
  double sum_tmp;
   
   X->sde->V[0](init, dest, X->sde->params);

  for (j=0; j<X->sde->dim_y; j++){

    sum_tmp=0.0;
    for (i=1; i<=X->sde->dim_BM; i++){

      X->sde->V[i](init, X->drift_step_interv, X->sde->params);
      X->sde->drift_corrector[i](init, X->drift_corrector_matrix, 
				 X->sde->params);      
      for (k=0; k<X->sde->dim_y; k++){
	sum_tmp +=0.5*X->drift_step_interv[k]*(X->drift_corrector_matrix[j][k]);
      } /* for (k) */
    } /* for (i) */
    dest[j]=dest[j]-sum_tmp;
   } /* for (j) */
}



int next_SDE_WA(SDE_WA_SLTN *X, double s, double y[], double dy[], void *rv){

  int j;
  double *sp_ptr;
  RV_NV *nv_sp_ptr;

  switch (X->alg){
  case 0:
    sp_ptr=(double *)rv;
    for (j=0; j<X->sde->dim_BM; j++) X->sample_pt.em[j]=sp_ptr[j];
    break;
  case 1:
    nv_sp_ptr=(RV_NV *)rv;
    X->sample_pt.nv->rv_nv_b=nv_sp_ptr->rv_nv_b;
    for (j=0; j<X->sde->dim_BM;j++)
      X->sample_pt.nv->rv_nv_n[j]=nv_sp_ptr->rv_nv_n[j];
    break;
  case 2:
    sp_ptr=(double *)rv;
    for (j=0; j<2*X->sde->dim_BM; j++) X->sample_pt.nn[j]=sp_ptr[j];
    break;
  }
  for (j=0; j<X->sde->dim_y; j++){
    X->initv[j]=y[j];
  }
  X->one_step(X, s);
  for (j=0; j<X->sde->dim_y; j++){
    dy[j]=X->destv[j];
  }
  return SDE_WA_SUCCESS;
}

