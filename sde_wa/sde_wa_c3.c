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
#include "sde_wa_butcher.h"
#include "sde_wa_c3.h"

/* exp(W_1)(x) */
int c3_vf_W(SDE_WA_SLTN *X, double s, double *step_initv, double *step_destv) 
{
	int i,j;
	double sq_s=sqrt(s);
  
	if (X->sde->sde_type==ITO){    	
		Ito_to_Strt_drift(X, step_initv, step_destv);
	}else {
		X->sde->V[0](step_initv, step_destv, X->sde->params);
	}

	for (j=0; j<X->sde->dim_y; j++){
		step_destv[j]=s*step_destv[j];
	}

	for (i=1; i<=X->sde->dim_BM; i++){
		X->sde->V[i](step_initv, X->drift_step_interv, X->sde->params);
		for (j=0; j<X->sde->dim_y; j++){
			step_destv[j] = step_destv[j] + sq_s*X->sample_pt.c3[i-1] * X->drift_step_interv[j];
		}/* for j */
	}/* for i */
  
  return SDE_WA_SUCCESS;
}

int sde_wa_c3(SDE_WA_SLTN *X, double s)
{
	//buthcer
	sde_wa_butcher(X, s, 1);

	return SDE_WA_SUCCESS;
}

