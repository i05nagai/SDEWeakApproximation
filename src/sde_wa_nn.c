/*
 * sde_wa_c3.c
 */

/*
 * $Author$
 * $Revision$
 * $Source$
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

#include <math.h>
#include "sde_wa.h"
#include "sde_wa_nn.h"
#include "sde_wa_butcher.h"

/*
 * X->sample_pt[2d]={eta^1_1, eta^2_1,..., eta^d_1, 
 *                   eta^1_2, eta^2_2, ..., eta^d_2}
 *
 * std_distr_to_Gaussian(X) makes 
 *   X->sample_pt[2d]={S^1_1,...S^d_1, 
 *     S^1_2,...,S^d_2}
 */
int std_distr_to_Gaussian(SDE_WA_SLTN *X){
  int i;
  double sq_2=sqrt(2.0);

  for (i=0; i<X->sde->dim_BM; i++){
    X->nn_sample_pt_interv[i]
      =X->sample_pt.nn[i]/2.0+X->sample_pt.nn[i+X->sde->dim_BM]/sq_2;
    X->nn_sample_pt_interv[i+X->sde->dim_BM]
      =X->sample_pt.nn[i]/2.0-X->sample_pt.nn[i+X->sde->dim_BM]/sq_2;
    
    X->sample_pt.nn[i]=X->nn_sample_pt_interv[i];
    X->sample_pt.nn[i+X->sde->dim_BM]=X->nn_sample_pt_interv[i+X->sde->dim_BM];
    
  }
  return SDE_WA_SUCCESS;
}

/* jth = 2 or 1 because exp(W_1)exp(W_2)(x) */
int nn_vf_W(SDE_WA_SLTN *X, double s, 
	    double *step_initv, double *step_destv, int jth){
  int i,j;
  double sq_s=sqrt(s);
  
  if (X->sde->sde_type==ITO){    	
    Ito_to_Strt_drift(X, step_initv, step_destv);
  }else {
    X->sde->V[0](step_initv, step_destv, X->sde->params);
  }

  for (j=0; j<X->sde->dim_y; j++){
    step_destv[j]=0.5*s*step_destv[j];
  }

  for (i=1; i<=X->sde->dim_BM; i++){
    X->sde->V[i](step_initv, X->drift_step_interv, X->sde->params);
    for (j=0; j<X->sde->dim_y; j++){
      step_destv[j]=step_destv[j]
	+sq_s*X->sample_pt.nn[(jth-1)*(X->sde->dim_BM)+i-1]
	*X->drift_step_interv[j];
    }/* for j */
  }/* for i */
  
  return SDE_WA_SUCCESS;
}

int sde_wa_nn(SDE_WA_SLTN *X, double s){
  double *tmp_nn_interv;
  
  /** std -> Gaussian **/
  std_distr_to_Gaussian(X);
  //sde_wa_butcher(x, s, 2);
  //buthcer or yoshida
  X->exp_eval(X, s, 2);
  
  tmp_nn_interv=X->initv;
  X->initv=X->destv;
  X->destv=tmp_nn_interv;
  //buthcer or yoshida
  X->exp_eval(X, s, 1);
  //sde_wa_butcher(X, s, 1);
  
  return SDE_WA_SUCCESS;
}
