/*
 * sde_wa_yoshida.c
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
#include <math.h>
#include "sde_wa_yoshida.h"
#include "sde_wa_butcher.h"

#define SOLUTION5B

#ifdef SOLUTION5A
#define W51 -0.117767998417887e1
#define W52 0.235573213359357
#define W53 0.784513610477560
#define W50 (1.0-2.0*(W51+W52+W53))
#endif

#ifdef SOLUTION5B
#define W51 -0.213228522200144e1
#define W52 0.426068187079180e-2
#define W53 0.143984816797668e1
#define W50 (1.0-2.0*(W51+W52+W53))
#endif

#ifdef SOLUTION5C
#define W51 0.152886228424922e-2
#define W52 -0.214403531630539e1
#define W53 0.144778256239930e1
#define W50 (1.0-2.0*(W51+W52+W53))
#endif

#define SOLUTION7A

#ifdef SOLUTION7A
#define W71 -0.161582374150097e1
#define W72 -0.244699182370524e1
#define W73 -0.716989419708120e-2
#define W74 0.244002732616735e1
#define W75 0.157739928123617
#define W76 0.182020630970714e1
#define W77 0.104242620869991e1
#define W70 (1.0-2.0*(W71+W72+W73+W74+W75+W76+W77))
#endif

#ifdef SOLUTION7B
#define W71 -0.169248587770116e-2
#define W72 0.289195744315849e1
#define W73 0.378039588360192e-2
#define W74 -0.289688250328827e1
#define W75 0.289105148970595e1
#define W76 -0.233864815101035
#define W77 0.148819229202922
#define W70 (1.0-2.0*(W71+W72+W73+W74+W75+W76+W77))
#endif

#ifdef SOLUTION7C
#define W71 0.311790812418427
#define W72 -0.155946803821447e1
#define W73 -0.167896928259640e1
#define W74 0.166335809963315e1
#define W75 -0.106458714789183e1
#define W76 0.1369349464169871e1
#define W77 0.629030650210433
#define W70 (1.0-2.0*(W71+W72+W73+W74+W75+W76+W77))
#endif

#ifdef SOLUTION7D
#define W71 0.102799849391985
#define W72 -0.196061023297549e1
#define W73 0.193813913762276e1
#define W74 -0.158240635368243
#define W75 -0.144485223686048e1
#define W76 0.253693336566229
#define W77 0.914844246229740
#define W70 (1.0-2.0*(W71+W72+W73+W74+W75+W76+W77))
#endif

#ifdef SOLUTION7E
#define W71 0.227738840094906e-1
#define W72 0.252778927322839e1
#define W73 -0.719180053552772e-1
#define W74 0.536018921307285e-2
#define W75 -0.204809795887393e1
#define W76 0.107990467703699
#define W77 0.130300165760014e1
#define W70 (1.0-2.0*(W71+W72+W73+W74+W75+W76+W77))
#endif
void sde_wa_yoshida(SDE_WA_SLTN *X, double s, int jth){
	double *tmp_interv;
  
	//switch (X->mth_is){
	switch (5){
	case 4:
		if (X->alg==N_V){ 
			return;
		}else if (X->alg==N_N){ //4th order symplectic integrator
			double W[] = {1.0, 0.5};
			exp_2th(X, s, jth, 0, W[1], W, 1);

			exp_2th(X, s, jth, 0, W[0], W, 1);

			exp_2th(X, s, jth, 0, W[1], W, 1);

			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else if (X->alg==C_3) {
			double W[] = {1.0, 0.5};
			exp_2th_c3(X, s, jth, 0, W[1], W, 1);

			exp_2th_c3(X, s, jth, 0, W[0], W, 1);

			exp_2th_c3(X, s, jth, 0, W[1], W, 1);

			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;

		}
		break;

	case 5:
		if (X->alg==N_V){ 
			return;
		}else if (X->alg==N_N){ //6th order symplectic integrator
			double W[] = {W50, W51, W52, W53};
//			exp_2th_Heston(X, s, jth, 0, W[3], W, 3);
//
//			exp_2th_Heston(X, s, jth, 0, W[2], W, 3);
//
//			exp_2th_Heston(X, s, jth, 0, W[1], W, 3);
//
//			exp_2th_Heston(X, s, jth, 0, W[0], W, 3);
//
//			exp_2th_Heston(X, s, jth, 0, W[1], W, 3);
//
//			exp_2th_Heston(X, s, jth, 0, W[2], W, 3);
//
//			exp_2th_Heston(X, s, jth, 0, W[3], W, 3);

//			exp_2th(X, s, jth, 0, W[3], W, 3);
//
//			exp_2th(X, s, jth, 0, W[2], W, 3);
//
//			exp_2th(X, s, jth, 0, W[1], W, 3);
//
//			exp_2th(X, s, jth, 0, W[0], W, 3);
//
//			exp_2th(X, s, jth, 0, W[1], W, 3);
//
//			exp_2th(X, s, jth, 0, W[2], W, 3);
//
//			exp_2th(X, s, jth, 0, W[3], W, 3);

			exp_2th_Kai(X, s, jth, 0, W[3], W, 3);

			exp_2th_Kai(X, s, jth, 0, W[2], W, 3);

			exp_2th_Kai(X, s, jth, 0, W[1], W, 3);

			exp_2th_Kai(X, s, jth, 0, W[0], W, 3);

			exp_2th_Kai(X, s, jth, 0, W[1], W, 3);

			exp_2th_Kai(X, s, jth, 0, W[2], W, 3);

			exp_2th_Kai(X, s, jth, 0, W[3], W, 3);

			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else if (X->alg==C_3) {
			double W[] = {W50, W51, W52, W53};
			exp_2th_c3(X, s, jth, 0, W[3], W, 3);

			exp_2th_c3(X, s, jth, 0, W[2], W, 3);

			exp_2th_c3(X, s, jth, 0, W[1], W, 3);

			exp_2th_c3(X, s, jth, 0, W[0], W, 3);

			exp_2th_c3(X, s, jth, 0, W[1], W, 3);

			exp_2th_c3(X, s, jth, 0, W[2], W, 3);

			exp_2th_c3(X, s, jth, 0, W[3], W, 3);
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;

		}
		break;
	case 7:
		if (X->alg==N_V){ 
			return;
		}else if(X->alg==N_N){ /* 8th order symplectic integrator */
			double W[] = {W70, W71, W72, W73, W74, W75, W76, W77};
			exp_2th(X, s, jth, 0, W[7], W, 7);

			exp_2th(X, s, jth, 0, W[6], W, 7);

			exp_2th(X, s, jth, 0, W[5], W, 7);

			exp_2th(X, s, jth, 0, W[4], W, 7);

			exp_2th(X, s, jth, 0, W[3], W, 7);

			exp_2th(X, s, jth, 0, W[2], W, 7);

			exp_2th(X, s, jth, 0, W[1], W, 7);

			exp_2th(X, s, jth, 0, W[0], W, 7);

			exp_2th(X, s, jth, 0, W[1], W, 7);

			exp_2th(X, s, jth, 0, W[2], W, 7);

			exp_2th(X, s, jth, 0, W[3], W, 7);

			exp_2th(X, s, jth, 0, W[4], W, 7);

			exp_2th(X, s, jth, 0, W[5], W, 7);

			exp_2th(X, s, jth, 0, W[6], W, 7);

			exp_2th(X, s, jth, 0, W[7], W, 7);

			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		}else if(X->alg==C_3) {/* else */
			double W[] = {W70, W71, W72, W73, W74, W75, W76, W77};
			exp_2th_c3(X, s, jth, 0, W[7], W, 7);

			exp_2th_c3(X, s, jth, 0, W[6], W, 7);

			exp_2th_c3(X, s, jth, 0, W[5], W, 7);

			exp_2th_c3(X, s, jth, 0, W[4], W, 7);

			exp_2th_c3(X, s, jth, 0, W[3], W, 7);

			exp_2th_c3(X, s, jth, 0, W[2], W, 7);

			exp_2th_c3(X, s, jth, 0, W[1], W, 7);

			exp_2th_c3(X, s, jth, 0, W[0], W, 7);

			exp_2th_c3(X, s, jth, 0, W[1], W, 7);

			exp_2th_c3(X, s, jth, 0, W[2], W, 7);

			exp_2th_c3(X, s, jth, 0, W[3], W, 7);

			exp_2th_c3(X, s, jth, 0, W[4], W, 7);

			exp_2th_c3(X, s, jth, 0, W[5], W, 7);

			exp_2th_c3(X, s, jth, 0, W[6], W, 7);

			exp_2th_c3(X, s, jth, 0, W[7], W, 7);

			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		}
		break;
	default:
		fprintf(stderr, "\nmth_is= 5 or 7\n");
		exit(0);
	}/* switch (m) */
}

/**
 * recursive 2nd order symplectic integrator for NN VF.
 * @param X SDE WA_SLTN.
 * @param s step size.
 * @param jth exp(W_jth) where W_j=0.5*s*V_0+sqrt(s)*Z^i_j*V_i.
 * @param d index of the VF.
 * @param w current weight.
 * @param W weight array.
 * @param length the length of W.
 * @return void. 
 */
void exp_2th(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length) {
	int i;
	double *tmp_interv;

	//V0
	if (d==0) {
		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/4.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/4.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;

		if (1==X->sde->dim_BM) {
			double sq_s=sqrt(s);
			if (X->sde->exp_sV[1]==NULL){
				sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], d+1);
			}else{
				X->sde->exp_sV[1](w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], X->initv, X->destv, X->sde->params);
			}
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else {
			for (i=length; i>=1; i--) {
				exp_2th(X, s, jth, 1, W[i]*w, W, length);
			}
			exp_2th(X, s, jth, 1, W[0]*w, W, length);
			for (i=1; i<=length; i++) {
				exp_2th(X, s, jth, 1, W[i]*w, W, length);
			}
		}

		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/4.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/4.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;

		return;
	}

	//Vi
	if (d<X->sde->dim_BM) {
		double sq_s=sqrt(s);
		if (X->sde->exp_sV[d]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, d);   
		}else{
			X->sde->exp_sV[d](w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi1:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		if (d==X->sde->dim_BM-1) {
			double sq_s=sqrt(s);
			if (X->sde->exp_sV[d+1]==NULL){
				sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], d+1);   
			}else{
				X->sde->exp_sV[d+1](w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], X->initv, X->destv, X->sde->params);
			}
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else {
			for (i=length; i>=1; i--) {
				exp_2th(X, s, jth, d+1, W[i]*w, W, length);
			}
			exp_2th(X, s, jth, d+1, W[0]*w, W, length);
			for (i=1; i<=length; i++) {
				exp_2th(X, s, jth, d+1, W[i]*w, W, length);
			}
		}

		if (X->sde->exp_sV[d]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, d);   
		}else{
			X->sde->exp_sV[d](w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi2:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		return;
	}
}

void exp_2th_Omega(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length) {
	int i;
	double *tmp_interv;

	//V0
	if (d==0 || d==2) {
		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/4.0/3.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/4.0/3.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;

		for (i=length; i>=1; i--) {
			exp_2th_Omega(X, s, jth, d+1, W[i]*w, W, length);
		}
		exp_2th_Omega(X, s, jth, d+1, W[0]*w, W, length);
		for (i=1; i<=length; i++) {
			exp_2th_Omega(X, s, jth, d+1, W[i]*w, W, length);
		}

		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/4.0/3.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/4.0/3.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;

		return;
	}

	//Vi
	if (d==1 || d==3) {
		int bm = d;
		if (d==3) {
			bm=2;
		}
		double sq_s=sqrt(s);
		if (X->sde->exp_sV[bm]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(bm-1)+X->sde->dim_BM*(jth-1)]/2.0, bm);   
		}else{
			X->sde->exp_sV[bm](w*sq_s*X->sample_pt.nn[(bm-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi1:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		if (d==3) {
			if (X->sde->exp_sV[0]==NULL){
				sde_wa_butcher(X, w*s/4.0/3.0, 0);   
			}else{
				X->sde->exp_sV[0](w*s/4.0/3.0, X->initv, X->destv, X->sde->params);
			}
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else {
			for (i=length; i>=1; i--) {
				exp_2th_Omega(X, s, jth, d+1, W[i]*w, W, length);
			}
			exp_2th_Omega(X, s, jth, d+1, W[0]*w, W, length);
			for (i=1; i<=length; i++) {
				exp_2th_Omega(X, s, jth, d+1, W[i]*w, W, length);
			}
		}

		if (X->sde->exp_sV[bm]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(bm-1)+X->sde->dim_BM*(jth-1)]/2.0, bm);   
		}else{
			X->sde->exp_sV[bm](w*sq_s*X->sample_pt.nn[(bm-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi2:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		return;
	}
}

void exp_2th_Kai(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length) {
	int i;
	double *tmp_interv;

	//V0
	if (d==0) {
		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/4.0/2.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/4.0/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;

		if (1==X->sde->dim_BM) {
			double sq_s=sqrt(s);
			if (X->sde->exp_sV[1]==NULL){
				sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], d+1);
			}else{
				X->sde->exp_sV[1](w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], X->initv, X->destv, X->sde->params);
			}
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else {
			for (i=length; i>=1; i--) {
				exp_2th_Kai(X, s, jth, d+1, W[i]*w, W, length);
			}
			exp_2th_Kai(X, s, jth, d+1, W[0]*w, W, length);
			for (i=1; i<=length; i++) {
				exp_2th_Kai(X, s, jth, d+1, W[i]*w, W, length);
			}
		}

		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/4.0/2.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/4.0/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;

		return;
	}

	//Vi
	if (d<X->sde->dim_BM) {
		double sq_s=sqrt(s);
		if (X->sde->exp_sV[d]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, d);   
		}else{
			X->sde->exp_sV[d](w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi1:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		if (d==X->sde->dim_BM) {
			if (X->sde->exp_sV[0]==NULL){
				sde_wa_butcher(X, w*s/4.0/2.0, 0);   
			}else{
				X->sde->exp_sV[0](w*s/4.0/2.0, X->initv, X->destv, X->sde->params);
			}
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else {
			for (i=length; i>=1; i--) {
				exp_2th_Kai(X, s, jth, d+1, W[i]*w, W, length);
			}
			exp_2th_Kai(X, s, jth, d+1, W[0]*w, W, length);
			for (i=1; i<=length; i++) {
				exp_2th_Kai(X, s, jth, d+1, W[i]*w, W, length);
			}
		}

		if (X->sde->exp_sV[d]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, d);   
		}else{
			X->sde->exp_sV[d](w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi2:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		return;
	}
}

void exp_2th_Heston(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length) {
	int i;
	double *tmp_interv;

	//V0
	if (d==0) {
		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/4.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/4.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("V01:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		for (i=length; i>=1; i--) {
			exp_2th_Heston(X, s, jth, 1, W[i]*w, W, length);
		}
		exp_2th_Heston(X, s, jth, 1, W[0]*w, W, length);
		for (i=1; i<=length; i++) {
			exp_2th_Heston(X, s, jth, 1, W[i]*w, W, length);
		}

		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/4.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/4.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("V02:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		return;
	}

	//V0
	if (d==1) {
		X->sde->exp_sV[X->sde->dim_BM+1](w*s/4.0, X->initv, X->destv, X->sde->params);

		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("V01:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		for (i=length; i>=1; i--) {
			exp_2th_Heston(X, s, jth, 2, W[i]*w, W, length);
		}
		exp_2th_Heston(X, s, jth, 2, W[0]*w, W, length);
		for (i=1; i<=length; i++) {
			exp_2th_Heston(X, s, jth, 2, W[i]*w, W, length);
		}

		X->sde->exp_sV[X->sde->dim_BM+1](w*s/4.0, X->initv, X->destv, X->sde->params);

		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("V02:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		return;
	}

	//Vi
	if (d<=X->sde->dim_BM+1) {
		d=d-1;
		double sq_s=sqrt(s);
		if (X->sde->exp_sV[d]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, d);   
		}else{
			X->sde->exp_sV[d](w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi1:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		if (d-1==X->sde->dim_BM) {
			double sq_s=sqrt(s);
			if (X->sde->exp_sV[d+1]==NULL){
				sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], d+1);   
			}else{
				X->sde->exp_sV[d+1](w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], X->initv, X->destv, X->sde->params);
			}
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else {
			for (i=length; i>=1; i--) {
				exp_2th_Heston(X, s, jth, d+2, W[i]*w, W, length);
			}
			exp_2th_Heston(X, s, jth, d+2, W[0]*w, W, length);
			for (i=1; i<=length; i++) {
				exp_2th_Heston(X, s, jth, d+2, W[i]*w, W, length);
			}
		}

		if (X->sde->exp_sV[d]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, d);   
		}else{
			X->sde->exp_sV[d](w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi2:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		return;
	}
}


/**
 * recursive 2nd order symplectic integrator for C3 VF.
 * @param X SDE WA_SLTN.
 * @param s step size.
 * @param jth exp(W_jth) where W_j=0.5*s*V_0+sqrt(s)*Z^i_j*V_i.
 * @param d index of the VF.
 * @param w current weight.
 * @param W weight array.
 * @param length the length of W.
 * @return void. 
 */
void exp_2th_c3(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length) {
	int i;
	double *tmp_interv;

	//V0
	if (d==0) {
		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/2.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("V01:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		if (1==X->sde->dim_BM) {
			double sq_s=sqrt(s);
			if (X->sde->exp_sV[d+1]==NULL){
				sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], d+1);
			}else{
				X->sde->exp_sV[d+1](w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], X->initv, X->destv, X->sde->params);
			}
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else {
			for (i=length; i>=1; i--) {
				exp_2th(X, s, jth, 1, W[i]*w, W, length);
			}
			exp_2th(X, s, jth, 1, W[0]*w, W, length);
			for (i=1; i<=length; i++) {
				exp_2th(X, s, jth, 1, W[i]*w, W, length);
			}
		}

		if (X->sde->exp_sV[0]==NULL){
			sde_wa_butcher(X, w*s/2.0, 0);   
		}else{
			X->sde->exp_sV[0](w*s/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("V02:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		return;
	}

	//Vi
	if (d<X->sde->dim_BM) {
		double sq_s=sqrt(s);
		if (X->sde->exp_sV[d]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, d);   
		}else{
			X->sde->exp_sV[d](w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi1:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		if (d==X->sde->dim_BM-1) {
			double sq_s=sqrt(s);
			if (X->sde->exp_sV[d+1]==NULL){
				sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], d+1);   
			}else{
				X->sde->exp_sV[d+1](w*sq_s*X->sample_pt.nn[d+X->sde->dim_BM*(jth-1)], X->initv, X->destv, X->sde->params);
			}
			tmp_interv=X->initv;
			X->initv=X->destv;
			X->destv=tmp_interv;
		} else {
			for (i=length; i>=1; i--) {
				exp_2th(X, s, jth, d+1, W[i]*w, W, length);
			}
			exp_2th(X, s, jth, d+1, W[0]*w, W, length);
			for (i=1; i<=length; i++) {
				exp_2th(X, s, jth, d+1, W[i]*w, W, length);
			}
		}

		if (X->sde->exp_sV[d]==NULL){
			sde_wa_butcher(X, w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, d);   
		}else{
			X->sde->exp_sV[d](w*sq_s*X->sample_pt.nn[(d-1)+X->sde->dim_BM*(jth-1)]/2.0, X->initv, X->destv, X->sde->params);
		}
		tmp_interv=X->initv;
		X->initv=X->destv;
		X->destv=tmp_interv;
		//printf("Vi2:%lf %lf %lf\n", X->initv[0], X->initv[1], X->initv[2]);

		return;
	}
}
