/*
 * sde_wa_butcher.c
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

#include <stdlib.h>
#include "sde_wa_butcher.h"
#include "sde_wa_nn.h"
#include "sde_wa_c3.h"

void sde_wa_butcher(SDE_WA_SLTN *X, double s, int jth){
  
	double *Y1, *Y2, *Y3, *Y4, *Y5;
	double *fY0, *fY1, *fY2, *fY3, *fY4, *fY5;
	int pf;
	int i;
  
	switch (X->mth_is){
	case 5:
    Y1 = X->rk_step_interv;
    Y2 = Y1+X->sde->dim_y;
    Y3 = Y1+2*X->sde->dim_y;
    Y4 = Y1+3*X->sde->dim_y;
    Y5 = Y1+4*X->sde->dim_y;
    fY0 = Y1+5*X->sde->dim_y;
    fY1 = Y1+6*X->sde->dim_y;
    fY2 = Y1+7*X->sde->dim_y;
    fY3 = Y1+8*X->sde->dim_y;
    fY4 = Y1+9*X->sde->dim_y;
    fY5 = Y1+10*X->sde->dim_y;

	if (X->alg==N_V){ /* 5th RK for N_V */
		if (jth==0 && X->sde->sde_type==ITO){/* V_0 && Ito form */
			Ito_to_Strt_drift(X, X->initv, fY0);
			for (i=0; i<X->sde->dim_y; i++)
			Y1[i] = X->initv[i] + (2.0/5.0)*s*fY0[i]; 
			Ito_to_Strt_drift(X, Y1, fY1);
			for (i=0; i<X->sde->dim_y; i++)
			Y2[i] = X->initv[i] + (11.0/64.0)*s*fY0[i] + (5.0/64.0)*s*fY1[i];
			Ito_to_Strt_drift(X, Y2, fY2);
			for (i=0; i<X->sde->dim_y; i++)
			Y3[i] = X->initv[i] + (1.0/2.0)*s*fY2[i];
			Ito_to_Strt_drift(X, Y3, fY3);
			for (i=0; i<X->sde->dim_y; i++)
			Y4[i] = X->initv[i] + (3.0/64.0)*s*fY0[i] - (15.0/64.0)*s*fY1[i]
			+ (3.0/8.0)*s*fY2[i] + (9.0/16.0)*s*fY3[i];
			Ito_to_Strt_drift(X, Y4, fY4);
			for (i=0; i<X->sde->dim_y; i++)
			Y5[i] = X->initv[i] + (5.0/7.0)*s*fY1[i] + (6.0/7.0)*s*fY2[i]
			-(12.0/7.0)*s*fY3[i] + (8.0/7.0)*s*fY4[i];
			Ito_to_Strt_drift(X, Y5, fY5);
			for (i=0; i<X->sde->dim_y; i++)
			X->destv[i] = X->initv[i] + (7.0/90.0)*s*fY0[i] + (32.0/90.0)*s*fY2[i]
			+ (12.0/90.0)*s*fY3[i] + (32.0/90.0)*s*fY4[i] 
			+ (7.0/90.0)*s*fY5[i];	
		}else{ /* V_1,..., V_d or Stratonovich form */
			pf = X->sde->V[jth](X->initv, fY0, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y1[i] = X->initv[i] + (2.0/5.0)*s*fY0[i]; 
			pf = X->sde->V[jth](Y1, fY1, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y2[i] = X->initv[i] + (11.0/64.0)*s*fY0[i] + (5.0/64.0)*s*fY1[i];
			pf = X->sde->V[jth](Y2, fY2, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y3[i] = X->initv[i] + (1.0/2.0)*s*fY2[i];
			pf = X->sde->V[jth](Y3, fY3, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y4[i] = X->initv[i] + (3.0/64.0)*s*fY0[i] - (15.0/64.0)*s*fY1[i]
			+ (3.0/8.0)*s*fY2[i] + (9.0/16.0)*s*fY3[i];
			pf = X->sde->V[jth](Y4, fY4, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y5[i] = X->initv[i] + (5.0/7.0)*s*fY1[i] + (6.0/7.0)*s*fY2[i]
			-(12.0/7.0)*s*fY3[i] + (8.0/7.0)*s*fY4[i];
			pf = X->sde->V[jth](Y5, fY5, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			X->destv[i] = X->initv[i] + (7.0/90.0)*s*fY0[i] + (32.0/90.0)*s*fY2[i]
			+ (12.0/90.0)*s*fY3[i] + (32.0/90.0)*s*fY4[i] + (7.0/90.0)*s*fY5[i];
		}
	} else if(X->alg==N_N) { /* 5th RK for N_N */
		// 5th RK for yoshida algorithm
		if (X->exp_type==YOSHIDA) {
			if (jth==0 && X->sde->sde_type==ITO){/* V_0 && Ito form */
				Ito_to_Strt_drift(X, X->initv, fY0);
				for (i=0; i<X->sde->dim_y; i++) {
					Y1[i] = X->initv[i] + (2.0/5.0)*s*fY0[i]; 
				}
				Ito_to_Strt_drift(X, Y1, fY1);
				for (i=0; i<X->sde->dim_y; i++) {
					Y2[i] = X->initv[i] + (11.0/64.0)*s*fY0[i] + (5.0/64.0)*s*fY1[i];
				}
				Ito_to_Strt_drift(X, Y2, fY2);
				for (i=0; i<X->sde->dim_y; i++) {
					Y3[i] = X->initv[i] + (1.0/2.0)*s*fY2[i];
				}
				Ito_to_Strt_drift(X, Y3, fY3);
				for (i=0; i<X->sde->dim_y; i++) {
					Y4[i] = X->initv[i] + (3.0/64.0)*s*fY0[i] - (15.0/64.0)*s*fY1[i]
						+ (3.0/8.0)*s*fY2[i] + (9.0/16.0)*s*fY3[i];
				}
				Ito_to_Strt_drift(X, Y4, fY4);
				for (i=0; i<X->sde->dim_y; i++) {
					Y5[i] = X->initv[i] + (5.0/7.0)*s*fY1[i] + (6.0/7.0)*s*fY2[i]
						-(12.0/7.0)*s*fY3[i] + (8.0/7.0)*s*fY4[i];
				}
				Ito_to_Strt_drift(X, Y5, fY5);
				for (i=0; i<X->sde->dim_y; i++) {
					X->destv[i] = X->initv[i] + (7.0/90.0)*s*fY0[i] + (32.0/90.0)*s*fY2[i]
						+ (12.0/90.0)*s*fY3[i] + (32.0/90.0)*s*fY4[i] 
						+ (7.0/90.0)*s*fY5[i];	
				}
			}else{ /* V_1,..., V_d or Stratonovich form */
				pf = X->sde->V[jth](X->initv, fY0, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y1[i] = X->initv[i] + (2.0/5.0)*s*fY0[i]; 
				}
				pf = X->sde->V[jth](Y1, fY1, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y2[i] = X->initv[i] + (11.0/64.0)*s*fY0[i] + (5.0/64.0)*s*fY1[i];
				}
				pf = X->sde->V[jth](Y2, fY2, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y3[i] = X->initv[i] + (1.0/2.0)*s*fY2[i];
				}
				pf = X->sde->V[jth](Y3, fY3, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y4[i] = X->initv[i] + (3.0/64.0)*s*fY0[i] - (15.0/64.0)*s*fY1[i]
						+ (3.0/8.0)*s*fY2[i] + (9.0/16.0)*s*fY3[i];
				}
				pf = X->sde->V[jth](Y4, fY4, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y5[i] = X->initv[i] + (5.0/7.0)*s*fY1[i] + (6.0/7.0)*s*fY2[i]
						-(12.0/7.0)*s*fY3[i] + (8.0/7.0)*s*fY4[i];
				}
				pf = X->sde->V[jth](Y5, fY5, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					X->destv[i] = X->initv[i] + (7.0/90.0)*s*fY0[i] + (32.0/90.0)*s*fY2[i]
						+ (12.0/90.0)*s*fY3[i] + (32.0/90.0)*s*fY4[i] + (7.0/90.0)*s*fY5[i];
				}
			}
		} else if (X->alg==C_3) {
			pf = nn_vf_W(X, s, X->initv, fY0, jth);
			for (i=0; i<X->sde->dim_y; i++)
			Y1[i] = X->initv[i] + (2.0/5.0)*fY0[i]; 
			pf = nn_vf_W(X, s, Y1, fY1, jth);
			for (i=0; i<X->sde->dim_y; i++)
			Y2[i] = X->initv[i] + (11.0/64.0)*fY0[i] + (5.0/64.0)*fY1[i];
			pf = nn_vf_W(X, s, Y2,  fY2, jth);
			for (i=0; i<X->sde->dim_y; i++)
			Y3[i] = X->initv[i] + (1.0/2.0)*fY2[i];
			pf = nn_vf_W(X, s, Y3,  fY3, jth);
			for (i=0; i<X->sde->dim_y; i++)
			Y4[i] = X->initv[i] + (3.0/64.0)*fY0[i] - (15.0/64.0)*fY1[i]
			+ (3.0/8.0)*fY2[i] + (9.0/16.0)*fY3[i];
			pf = nn_vf_W(X, s, Y4,  fY4, jth);
			for (i=0; i<X->sde->dim_y; i++)
			Y5[i] = X->initv[i] + (5.0/7.0)*fY1[i] + (6.0/7.0)*fY2[i]
			-(12.0/7.0)*fY3[i] + (8.0/7.0)*fY4[i];
			pf = nn_vf_W(X, s, Y5,  fY5, jth);
			for (i=0; i<X->sde->dim_y; i++)
			X->destv[i] = X->initv[i] + (7.0/90.0)*fY0[i] + (32.0/90.0)*fY2[i]
			+ (12.0/90.0)*fY3[i] + (32.0/90.0)*fY4[i] + (7.0/90.0)*fY5[i];  
		}
    } else if (X->alg==C_3) { /* 5th RK for Cubature 3 */
		// 5th RK for yoshida algorithm
		if (X->exp_type==YOSHIDA) {
			if (jth==0 && X->sde->sde_type==ITO){/* V_0 && Ito form */
				Ito_to_Strt_drift(X, X->initv, fY0);
				for (i=0; i<X->sde->dim_y; i++) {
					Y1[i] = X->initv[i] + (2.0/5.0)*s*fY0[i]; 
				}
				Ito_to_Strt_drift(X, Y1, fY1);
				for (i=0; i<X->sde->dim_y; i++) {
					Y2[i] = X->initv[i] + (11.0/64.0)*s*fY0[i] + (5.0/64.0)*s*fY1[i];
				}
				Ito_to_Strt_drift(X, Y2, fY2);
				for (i=0; i<X->sde->dim_y; i++) {
					Y3[i] = X->initv[i] + (1.0/2.0)*s*fY2[i];
				}
				Ito_to_Strt_drift(X, Y3, fY3);
				for (i=0; i<X->sde->dim_y; i++) {
					Y4[i] = X->initv[i] + (3.0/64.0)*s*fY0[i] - (15.0/64.0)*s*fY1[i]
						+ (3.0/8.0)*s*fY2[i] + (9.0/16.0)*s*fY3[i];
				}
				Ito_to_Strt_drift(X, Y4, fY4);
				for (i=0; i<X->sde->dim_y; i++) {
					Y5[i] = X->initv[i] + (5.0/7.0)*s*fY1[i] + (6.0/7.0)*s*fY2[i]
						-(12.0/7.0)*s*fY3[i] + (8.0/7.0)*s*fY4[i];
				}
				Ito_to_Strt_drift(X, Y5, fY5);
				for (i=0; i<X->sde->dim_y; i++) {
					X->destv[i] = X->initv[i] + (7.0/90.0)*s*fY0[i] + (32.0/90.0)*s*fY2[i]
						+ (12.0/90.0)*s*fY3[i] + (32.0/90.0)*s*fY4[i] 
						+ (7.0/90.0)*s*fY5[i];	
				}
			}else{ /* V_1,..., V_d or Stratonovich form */
				pf = X->sde->V[jth](X->initv, fY0, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y1[i] = X->initv[i] + (2.0/5.0)*s*fY0[i]; 
				}
				pf = X->sde->V[jth](Y1, fY1, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y2[i] = X->initv[i] + (11.0/64.0)*s*fY0[i] + (5.0/64.0)*s*fY1[i];
				}
				pf = X->sde->V[jth](Y2, fY2, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y3[i] = X->initv[i] + (1.0/2.0)*s*fY2[i];
				}
				pf = X->sde->V[jth](Y3, fY3, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y4[i] = X->initv[i] + (3.0/64.0)*s*fY0[i] - (15.0/64.0)*s*fY1[i]
						+ (3.0/8.0)*s*fY2[i] + (9.0/16.0)*s*fY3[i];
				}
				pf = X->sde->V[jth](Y4, fY4, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y5[i] = X->initv[i] + (5.0/7.0)*s*fY1[i] + (6.0/7.0)*s*fY2[i]
						-(12.0/7.0)*s*fY3[i] + (8.0/7.0)*s*fY4[i];
				}
				pf = X->sde->V[jth](Y5, fY5, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					X->destv[i] = X->initv[i] + (7.0/90.0)*s*fY0[i] + (32.0/90.0)*s*fY2[i]
						+ (12.0/90.0)*s*fY3[i] + (32.0/90.0)*s*fY4[i] + (7.0/90.0)*s*fY5[i];
				}
			}
		} else {
			pf = c3_vf_W(X, s, X->initv, fY0);
			for (i=0; i<X->sde->dim_y; i++) {
				Y1[i] = X->initv[i] + (2.0/5.0)*fY0[i]; 
			}
			pf = c3_vf_W(X, s, Y1, fY1);
			for (i=0; i<X->sde->dim_y; i++) {
				Y2[i] = X->initv[i] + (11.0/64.0)*fY0[i] + (5.0/64.0)*fY1[i];
			}
			pf = c3_vf_W(X, s, Y2,  fY2);
			for (i=0; i<X->sde->dim_y; i++) {
				Y3[i] = X->initv[i] + (1.0/2.0)*fY2[i];
			}
			pf = c3_vf_W(X, s, Y3,  fY3);
			for (i=0; i<X->sde->dim_y; i++) {
				Y4[i] = X->initv[i] + (3.0/64.0)*fY0[i] - (15.0/64.0)*fY1[i]
					+ (3.0/8.0)*fY2[i] + (9.0/16.0)*fY3[i];
			}
			pf = c3_vf_W(X, s, Y4,  fY4);
			for (i=0; i<X->sde->dim_y; i++) {
				Y5[i] = X->initv[i] + (5.0/7.0)*fY1[i] + (6.0/7.0)*fY2[i]
					-(12.0/7.0)*fY3[i] + (8.0/7.0)*fY4[i];
			}
			pf = c3_vf_W(X, s, Y5,  fY5);
			for (i=0; i<X->sde->dim_y; i++) {
				X->destv[i] = X->initv[i] + (7.0/90.0)*fY0[i] + (32.0/90.0)*fY2[i]
					+ (12.0/90.0)*fY3[i] + (32.0/90.0)*fY4[i] + (7.0/90.0)*fY5[i];  
			}

		}

	}
    break;
	case 7:
    {
	const double a21 = 1.0/6.0;
	const double a32 = 1.0/3.0;
	const double a41 = 1.0/8.0;
	const double a43 = 3.0/8.0;
	const double a51 = 148.0/1331.0;
	const double a53 = 150.0/1331.0;
	const double a54 = -56.0/1331.0;
	const double a61 = -404.0/243.0;
	const double a63 = -170.0/27.0;
	const double a64 = 4024.0/1701.0;
	const double a65 = 10648.0/1701.0;
	const double a71 = 2466.0/2401.0;
	const double a73 = 1242.0/343.0;
	const double a74 = -19176.0/16807.0;
	const double a75 = -51909.0/16807.0;
	const double a76 = 1053.0/2401.0;
	const double a81 = 5.0/154.0;
	const double a84 = 96.0/539.0;
	const double a85 = -1815.0/20384.0;
	const double a86 = -405.0/2464.0;
	const double a87 = 49.0/1144.0;
	const double a91 = -113.0/32.0;
	const double a93 = -195.0/22.0;
	const double a94 = 32.0/7.0;
	const double a95 = 29403.0/3584.0;
	const double a96 = -729.0/512.0;
	const double a97 = 1029.0/1408.0;
	const double a98 = 21.0/16.0;
	const double b4 = 32.0/105.0;
	const double b5 = 1771561.0/6289920.0;
	const double b6 = 243.0/2560.0;
	const double b7 = 16807.0/74880.0;
	const double b8 = 77.0/1440.0;
	const double b9 = 11.0/270.0;

	double *Y6, *Y7, *Y8;
	double *fY6, *fY7, *fY8;

	Y1=X->rk_step_interv;
	Y2=Y1+X->sde->dim_y;
	Y3=Y1+2*X->sde->dim_y;
	Y4=Y1+3*X->sde->dim_y;
	Y5=Y1+4*X->sde->dim_y;
	Y6=Y1+5*X->sde->dim_y;
	Y7=Y1+6*X->sde->dim_y;
	Y8=Y1+7*X->sde->dim_y;
	fY0=Y1+8*X->sde->dim_y;
	fY1=Y1+9*X->sde->dim_y;
	fY2=Y1+10*X->sde->dim_y;
	fY3=Y1+11*X->sde->dim_y;
	fY4=Y1+12*X->sde->dim_y;
	fY5=Y1+13*X->sde->dim_y;
	fY6=Y1+14*X->sde->dim_y;
	fY7=Y1+15*X->sde->dim_y;
	fY8=Y1+16*X->sde->dim_y;
      	  
	if (X->alg==N_V){ /* 7th RK for N_V */
		if (jth==0 && X->sde->sde_type==ITO){/* V_0 && Ito form */
			Ito_to_Strt_drift(X, X->initv, fY0);
			for (i=0; i<X->sde->dim_y; i++)
			Y1[i] = X->initv[i] + a21*s*fY0[i];
			Ito_to_Strt_drift(X, Y1, fY1);
			for (i=0; i<X->sde->dim_y; i++)
			Y2[i] = X->initv[i] + a32*s*fY1[i];
			Ito_to_Strt_drift(X, Y2, fY2);
			for (i=0; i<X->sde->dim_y; i++)
			Y3[i] = X->initv[i] + a41*s*fY0[i]+a43*s*fY2[i];
			Ito_to_Strt_drift(X, Y3, fY3);
			for (i=0; i<X->sde->dim_y; i++)
			Y4[i] = X->initv[i] + a51*s*fY0[i] + a53*s*fY2[i] + a54*s*fY3[i];
			Ito_to_Strt_drift(X, Y4, fY4);
			for (i=0; i<X->sde->dim_y; i++)
			Y5[i] = X->initv[i] + a61*s*fY0[i]
			  + a63*s*fY2[i] + a64*s*fY3[i] + a65*s*fY4[i];
			Ito_to_Strt_drift(X, Y5, fY5);
			for (i=0; i<X->sde->dim_y; i++)
			Y6[i] = X->initv[i] + a71*s*fY0[i]
			  + a73*s*fY2[i] + a74*s*fY3[i] + a75*s*fY4[i] + a76*s*fY5[i];
			Ito_to_Strt_drift(X, Y6, fY6);
			for (i=0; i<X->sde->dim_y; i++)
			Y7[i] = X->initv[i] + a81*s*fY0[i] + a84*s*fY3[i]
			  + a85*s*fY4[i] + a86*s*fY5[i] + a87*s*fY6[i];
			Ito_to_Strt_drift(X, Y7, fY7);
			for (i=0; i<X->sde->dim_y; i++)
			Y8[i] = X->initv[i] + a91*s*fY0[i] + +a93*s*fY2[i] + a94*s*fY3[i] 
			  + a95*s*fY4[i] 
			  + a96*s*fY5[i] +a97*s*fY6[i] + a98*s*fY7[i];
			Ito_to_Strt_drift(X, Y8, fY8);
			for (i=0; i<X->sde->dim_y; i++)
			X->destv[i] = X->initv[i] + b4*s*fY3[i] + b5*s*fY4[i] + b6*s*fY5[i] 
			  + b7*s*fY6[i] + b8*s*fY7[i] + b9*s*fY8[i];
		}else{ /* V_1, ... V_d, or Stratonovich form */
			pf=X->sde->V[jth](X->initv,fY0, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y1[i] = X->initv[i] + a21*s*fY0[i]; 
			pf = X->sde->V[jth](Y1, fY1, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y2[i] = X->initv[i] + a32*s*fY1[i];
			pf = X->sde->V[jth](Y2, fY2, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y3[i] = X->initv[i] + a41*s*fY0[i]+a43*s*fY2[i];
			pf = X->sde->V[jth](Y3, fY3, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y4[i] = X->initv[i] + a51*s*fY0[i] + a53*s*fY2[i] + a54*s*fY3[i];
			pf = X->sde->V[jth](Y4, fY4, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y5[i] = X->initv[i] + a61*s*fY0[i] + a63*s*fY2[i] + a64*s*fY3[i]
			+ a65*s*fY4[i];
			pf = X->sde->V[jth](Y5, fY5, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y6[i] = X->initv[i] + a71*s*fY0[i] + a73*s*fY2[i] + a74*s*fY3[i]
			+ a75*s*fY4[i] + a76*s*fY5[i];
			pf = X->sde->V[jth](Y6, fY6, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y7[i] = X->initv[i] + a81*s*fY0[i] + a84*s*fY3[i] + a85*s*fY4[i]
			+ a86*s*fY5[i] + a87*s*fY6[i];
			pf = X->sde->V[jth](Y7, fY7, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			Y8[i] = X->initv[i] + a91*s*fY0[i] + +a93*s*fY2[i] + a94*s*fY3[i] 
			+ a95*s*fY4[i] + a96*s*fY5[i] +a97*s*fY6[i] + a98*s*fY7[i];
			pf = X->sde->V[jth](Y8, fY8, X->sde->params);
			for (i=0; i<X->sde->dim_y; i++)
			X->destv[i] = X->initv[i] + b4*s*fY3[i] + b5*s*fY4[i] + b6*s*fY5[i] 
			+ b7*s*fY6[i] + b8*s*fY7[i] + b9*s*fY8[i];
		}
	}else if (X->alg==N_N) { /* 7th RK for N_N */
		// 7th RK for yoshida algorithm
		if (X->exp_type==YOSHIDA) {
			if (jth==0 && X->sde->sde_type==ITO){/* V_0 && Ito form */
				Ito_to_Strt_drift(X, X->initv, fY0);
				for (i=0; i<X->sde->dim_y; i++)
				Y1[i] = X->initv[i] + a21*s*fY0[i];
				Ito_to_Strt_drift(X, Y1, fY1);
				for (i=0; i<X->sde->dim_y; i++)
				Y2[i] = X->initv[i] + a32*s*fY1[i];
				Ito_to_Strt_drift(X, Y2, fY2);
				for (i=0; i<X->sde->dim_y; i++)
				Y3[i] = X->initv[i] + a41*s*fY0[i]+a43*s*fY2[i];
				Ito_to_Strt_drift(X, Y3, fY3);
				for (i=0; i<X->sde->dim_y; i++)
				Y4[i] = X->initv[i] + a51*s*fY0[i] + a53*s*fY2[i] + a54*s*fY3[i];
				Ito_to_Strt_drift(X, Y4, fY4);
				for (i=0; i<X->sde->dim_y; i++)
				Y5[i] = X->initv[i] + a61*s*fY0[i]
				  + a63*s*fY2[i] + a64*s*fY3[i] + a65*s*fY4[i];
				Ito_to_Strt_drift(X, Y5, fY5);
				for (i=0; i<X->sde->dim_y; i++)
				Y6[i] = X->initv[i] + a71*s*fY0[i]
				  + a73*s*fY2[i] + a74*s*fY3[i] + a75*s*fY4[i] + a76*s*fY5[i];
				Ito_to_Strt_drift(X, Y6, fY6);
				for (i=0; i<X->sde->dim_y; i++)
				Y7[i] = X->initv[i] + a81*s*fY0[i] + a84*s*fY3[i]
				  + a85*s*fY4[i] + a86*s*fY5[i] + a87*s*fY6[i];
				Ito_to_Strt_drift(X, Y7, fY7);
				for (i=0; i<X->sde->dim_y; i++)
				Y8[i] = X->initv[i] + a91*s*fY0[i] + +a93*s*fY2[i] + a94*s*fY3[i] 
				  + a95*s*fY4[i] 
				  + a96*s*fY5[i] +a97*s*fY6[i] + a98*s*fY7[i];
				Ito_to_Strt_drift(X, Y8, fY8);
				for (i=0; i<X->sde->dim_y; i++)
				X->destv[i] = X->initv[i] + b4*s*fY3[i] + b5*s*fY4[i] + b6*s*fY5[i] 
				  + b7*s*fY6[i] + b8*s*fY7[i] + b9*s*fY8[i];
			}else{ /* V_1, ... V_d, or Stratonovich form */
				pf=X->sde->V[jth](X->initv,fY0, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				Y1[i] = X->initv[i] + a21*s*fY0[i]; 
				pf = X->sde->V[jth](Y1, fY1, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				Y2[i] = X->initv[i] + a32*s*fY1[i];
				pf = X->sde->V[jth](Y2, fY2, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				Y3[i] = X->initv[i] + a41*s*fY0[i]+a43*s*fY2[i];
				pf = X->sde->V[jth](Y3, fY3, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				Y4[i] = X->initv[i] + a51*s*fY0[i] + a53*s*fY2[i] + a54*s*fY3[i];
				pf = X->sde->V[jth](Y4, fY4, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				Y5[i] = X->initv[i] + a61*s*fY0[i] + a63*s*fY2[i] + a64*s*fY3[i]
				+ a65*s*fY4[i];
				pf = X->sde->V[jth](Y5, fY5, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				Y6[i] = X->initv[i] + a71*s*fY0[i] + a73*s*fY2[i] + a74*s*fY3[i]
				+ a75*s*fY4[i] + a76*s*fY5[i];
				pf = X->sde->V[jth](Y6, fY6, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				Y7[i] = X->initv[i] + a81*s*fY0[i] + a84*s*fY3[i] + a85*s*fY4[i]
				+ a86*s*fY5[i] + a87*s*fY6[i];
				pf = X->sde->V[jth](Y7, fY7, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				Y8[i] = X->initv[i] + a91*s*fY0[i] + +a93*s*fY2[i] + a94*s*fY3[i] 
				+ a95*s*fY4[i] + a96*s*fY5[i] +a97*s*fY6[i] + a98*s*fY7[i];
				pf = X->sde->V[jth](Y8, fY8, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++)
				X->destv[i] = X->initv[i] + b4*s*fY3[i] + b5*s*fY4[i] + b6*s*fY5[i] 
				+ b7*s*fY6[i] + b8*s*fY7[i] + b9*s*fY8[i];
			}

		} else{ 
			pf=nn_vf_W(X, s, X->initv, fY0, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  Y1[i] = X->initv[i] + a21*fY0[i];
			pf = nn_vf_W(X, s, Y1, fY1, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  Y2[i] = X->initv[i] + a32*fY1[i];
			pf = nn_vf_W(X, s, Y2, fY2, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  Y3[i] = X->initv[i] + a41*fY0[i]+a43*fY2[i];
			pf = nn_vf_W(X, s, Y3, fY3, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  Y4[i] = X->initv[i] + a51*fY0[i] + a53*fY2[i] + a54*fY3[i];
			pf = nn_vf_W(X, s, Y4, fY4, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  Y5[i] = X->initv[i] + a61*fY0[i] + a63*fY2[i] + a64*fY3[i]
				+ a65*fY4[i];
			pf = nn_vf_W(X, s, Y5, fY5, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  Y6[i] = X->initv[i] + a71*fY0[i] + a73*fY2[i] + a74*fY3[i]
				+ a75*fY4[i] + a76*fY5[i];
			pf = nn_vf_W(X, s, Y6, fY6, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  Y7[i] = X->initv[i] + a81*fY0[i] + a84*fY3[i] + a85*fY4[i]
				+ a86*fY5[i] + a87*fY6[i];
			pf = nn_vf_W(X, s, Y7, fY7, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  Y8[i] = X->initv[i] + a91*fY0[i] + +a93*fY2[i] + a94*fY3[i] 
				+ a95*fY4[i] + a96*fY5[i] +a97*fY6[i] + a98*fY7[i];
			pf = nn_vf_W(X, s, Y8, fY8, jth);
			for (i=0; i<X->sde->dim_y; i++)
			  X->destv[i] = X->initv[i] + b4*fY3[i] + b5*fY4[i] + b6*fY5[i] 
				+ b7*fY6[i] + b8*fY7[i] + b9*fY8[i];
		}
	} else if (X->alg==C_3) { /* else */
		// 7th RK for yoshida algorithm
		if (X->exp_type==YOSHIDA) {
			if (jth==0 && X->sde->sde_type==ITO){/* V_0 && Ito form */
				Ito_to_Strt_drift(X, X->initv, fY0);
				for (i=0; i<X->sde->dim_y; i++)
				Y1[i] = X->initv[i] + a21*s*fY0[i];
				Ito_to_Strt_drift(X, Y1, fY1);
				for (i=0; i<X->sde->dim_y; i++)
				Y2[i] = X->initv[i] + a32*s*fY1[i];
				Ito_to_Strt_drift(X, Y2, fY2);
				for (i=0; i<X->sde->dim_y; i++)
				Y3[i] = X->initv[i] + a41*s*fY0[i]+a43*s*fY2[i];
				Ito_to_Strt_drift(X, Y3, fY3);
				for (i=0; i<X->sde->dim_y; i++)
				Y4[i] = X->initv[i] + a51*s*fY0[i] + a53*s*fY2[i] + a54*s*fY3[i];
				Ito_to_Strt_drift(X, Y4, fY4);
				for (i=0; i<X->sde->dim_y; i++)
				Y5[i] = X->initv[i] + a61*s*fY0[i]
				  + a63*s*fY2[i] + a64*s*fY3[i] + a65*s*fY4[i];
				Ito_to_Strt_drift(X, Y5, fY5);
				for (i=0; i<X->sde->dim_y; i++)
				Y6[i] = X->initv[i] + a71*s*fY0[i]
				  + a73*s*fY2[i] + a74*s*fY3[i] + a75*s*fY4[i] + a76*s*fY5[i];
				Ito_to_Strt_drift(X, Y6, fY6);
				for (i=0; i<X->sde->dim_y; i++)
				Y7[i] = X->initv[i] + a81*s*fY0[i] + a84*s*fY3[i]
				  + a85*s*fY4[i] + a86*s*fY5[i] + a87*s*fY6[i];
				Ito_to_Strt_drift(X, Y7, fY7);
				for (i=0; i<X->sde->dim_y; i++)
				Y8[i] = X->initv[i] + a91*s*fY0[i] + +a93*s*fY2[i] + a94*s*fY3[i] 
				  + a95*s*fY4[i] 
				  + a96*s*fY5[i] +a97*s*fY6[i] + a98*s*fY7[i];
				Ito_to_Strt_drift(X, Y8, fY8);
				for (i=0; i<X->sde->dim_y; i++)
				X->destv[i] = X->initv[i] + b4*s*fY3[i] + b5*s*fY4[i] + b6*s*fY5[i] 
				  + b7*s*fY6[i] + b8*s*fY7[i] + b9*s*fY8[i];
			}else{ /* V_1, ... V_d, or Stratonovich form */
				pf=X->sde->V[jth](X->initv,fY0, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y1[i] = X->initv[i] + a21*s*fY0[i]; 
				}
				pf = X->sde->V[jth](Y1, fY1, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y2[i] = X->initv[i] + a32*s*fY1[i];
				}
				pf = X->sde->V[jth](Y2, fY2, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y3[i] = X->initv[i] + a41*s*fY0[i]+a43*s*fY2[i];
				}
				pf = X->sde->V[jth](Y3, fY3, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y4[i] = X->initv[i] + a51*s*fY0[i] + a53*s*fY2[i] + a54*s*fY3[i];
				}
				pf = X->sde->V[jth](Y4, fY4, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y5[i] = X->initv[i] + a61*s*fY0[i] + a63*s*fY2[i] + a64*s*fY3[i]
						+ a65*s*fY4[i];
				}
				pf = X->sde->V[jth](Y5, fY5, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y6[i] = X->initv[i] + a71*s*fY0[i] + a73*s*fY2[i] + a74*s*fY3[i]
						+ a75*s*fY4[i] + a76*s*fY5[i];
				}
				pf = X->sde->V[jth](Y6, fY6, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y7[i] = X->initv[i] + a81*s*fY0[i] + a84*s*fY3[i] + a85*s*fY4[i]
						+ a86*s*fY5[i] + a87*s*fY6[i];
				}
				pf = X->sde->V[jth](Y7, fY7, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					Y8[i] = X->initv[i] + a91*s*fY0[i] + +a93*s*fY2[i] + a94*s*fY3[i] 
						+ a95*s*fY4[i] + a96*s*fY5[i] +a97*s*fY6[i] + a98*s*fY7[i];
				}
				pf = X->sde->V[jth](Y8, fY8, X->sde->params);
				for (i=0; i<X->sde->dim_y; i++) {
					X->destv[i] = X->initv[i] + b4*s*fY3[i] + b5*s*fY4[i] + b6*s*fY5[i] 
						+ b7*s*fY6[i] + b8*s*fY7[i] + b9*s*fY8[i];
				}
			}
		} else{ 
			pf = c3_vf_W(X, s, X->initv, fY0);
			for (i=0; i<X->sde->dim_y; i++) {
				Y1[i] = X->initv[i] + a21*fY0[i];
			}
			pf = c3_vf_W(X, s, Y1, fY1);
			for (i=0; i<X->sde->dim_y; i++) {
				Y2[i] = X->initv[i] + a32*fY1[i];
			}
			pf = c3_vf_W(X, s, Y2, fY2);
			for (i=0; i<X->sde->dim_y; i++) {
				Y3[i] = X->initv[i] + a41*fY0[i]+a43*fY2[i];
			}
			pf = c3_vf_W(X, s, Y3, fY3);
			for (i=0; i<X->sde->dim_y; i++) {
				Y4[i] = X->initv[i] + a51*fY0[i] + a53*fY2[i] + a54*fY3[i];
			}
			pf = c3_vf_W(X, s, Y4, fY4);
			for (i=0; i<X->sde->dim_y; i++) {
				Y5[i] = X->initv[i] + a61*fY0[i] + a63*fY2[i] + a64*fY3[i] + a65*fY4[i];
			}
			pf = c3_vf_W(X, s, Y5, fY5);
			for (i=0; i<X->sde->dim_y; i++) {
				Y6[i] = X->initv[i] + a71*fY0[i] + a73*fY2[i] + a74*fY3[i]
					+ a75*fY4[i] + a76*fY5[i];
			}
			pf = c3_vf_W(X, s, Y6, fY6);
			for (i=0; i<X->sde->dim_y; i++) {
				Y7[i] = X->initv[i] + a81*fY0[i] + a84*fY3[i] + a85*fY4[i]
					+ a86*fY5[i] + a87*fY6[i];
			}
			pf = c3_vf_W(X, s, Y7, fY7);
			for (i=0; i<X->sde->dim_y; i++) {
				Y8[i] = X->initv[i] + a91*fY0[i] + +a93*fY2[i] + a94*fY3[i] 
					+ a95*fY4[i] + a96*fY5[i] +a97*fY6[i] + a98*fY7[i];
			}
			pf = c3_vf_W(X, s, Y8, fY8);
			for (i=0; i<X->sde->dim_y; i++) {
				X->destv[i] = X->initv[i] + b4*fY3[i] + b5*fY4[i] + b6*fY5[i] 
					+ b7*fY6[i] + b8*fY7[i] + b9*fY8[i];
			}
		}

	}
	} /* const double ... */
    break;
  default:
    fprintf(stderr, "\nmth_is= 5 or 7\n");
    exit(0);
  }/* switch (m) */
}

