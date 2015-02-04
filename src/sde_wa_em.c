/*
 * sde_wa_em.c
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
#include "sde_wa_em.h"

void Strt_to_Ito_drift(SDE_WA_SLTN *X){

  int i,j,k;
  double sum_tmp;
   
  X->sde->V[0](X->initv, X->destv, X->sde->params);
  
  for (j=0; j<X->sde->dim_y; j++){
    sum_tmp=0.0;
    for (i=1; i<=X->sde->dim_BM; i++){
      X->sde->V[i](X->initv, X->drift_step_interv, X->sde->params);
      X->sde->drift_corrector[i](X->initv, X->drift_corrector_matrix, 
				       X->sde->params);
      for (k=0; k<X->sde->dim_y; k++){
	sum_tmp +=0.5*X->drift_step_interv[k]*(X->drift_corrector_matrix[j][k]);
      } /* for (k) */
    } /* for (i) */
    X->destv[j]=X->destv[j]+sum_tmp;
  } /* for (j) */
}

int sde_wa_em(SDE_WA_SLTN *X, double s){
  
  int j,i;
  double sq_s=sqrt(s);


  if (X->sde->sde_type==ITO){/* if Ito type */
    X->sde->V[0](X->initv, X->destv, X->sde->params);
  }else{/* if Stratonovich type*/
    Strt_to_Ito_drift(X);
  }

  for (j=0; j<X->sde->dim_y; j++){
    X->destv[j]=X->initv[j]+s*X->destv[j];
  }
  
  for (i=1; i<=X->sde->dim_BM; i++){
    X->sde->V[i](X->initv, X->drift_step_interv, X->sde->params);
    for (j=0;j<X->sde->dim_y; j++){
      X->destv[j]=X->destv[j]+X->drift_step_interv[j]*sq_s*X->sample_pt.em[i-1];
    }
  }
  return SDE_WA_SUCCESS;
}
