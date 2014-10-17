/*
 * sde_wa_nv.c
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
#include "sde_wa_nv.h"
#include "sde_wa_butcher.h"

int sde_wa_nv(SDE_WA_SLTN *X, double s){
  int i;
  double sq_s=sqrt(s);
  double *tmp_nv_interv;

  if (X->sde->exp_sV[0]==NULL){
    sde_wa_butcher(X, s/2.0, 0);
  }else{
    X->sde->exp_sV[0](s/2.0, X->initv, X->destv, X->sde->params);
  }

  tmp_nv_interv=X->initv;
  X->initv=X->destv;
  X->destv=tmp_nv_interv;

  if (X->sample_pt.nv->rv_nv_b==H){ /* if HEAD */
    for (i=1; i<=X->sde->dim_BM; i++){
      if (X->sde->exp_sV[i]==NULL){
	sde_wa_butcher(X, sq_s*X->sample_pt.nv->rv_nv_n[i-1], i);
      }else{
	X->sde->exp_sV[i](sq_s*X->sample_pt.nv->rv_nv_n[i-1], 
			  X->initv, X->destv, X->sde->params);
      }
      tmp_nv_interv=X->initv;
      X->initv=X->destv;
      X->destv=tmp_nv_interv;
    } /* for (i) */
  }else{ /* if TAIL */
    for (i=X->sde->dim_BM; i>=1; i--){
      if (X->sde->exp_sV[i]==NULL){
	sde_wa_butcher(X, sq_s*X->sample_pt.nv->rv_nv_n[i-1], i);   
      }else{
	X->sde->exp_sV[i]
	  (sq_s*X->sample_pt.nv->rv_nv_n[i-1], 
	   X->initv, X->destv, X->sde->params);
      }
      tmp_nv_interv=X->initv;
      X->initv=X->destv;
      X->destv=tmp_nv_interv;
    }/* for (i) */
  }
  if (X->sde->exp_sV[0]==NULL){
    sde_wa_butcher(X, s/2.0, 0);
  }else{
    X->sde->exp_sV[0](s/2.0, X->initv, X->destv, X->sde->params);
  }  
  return SDE_WA_SUCCESS;
}

