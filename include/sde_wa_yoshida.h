/*
 * sde_wa_yoshida.h
 */

#include "sde_wa.h"

#ifndef __SDE_WA_YOSHIDA_H__
#define __SDE_WA_YOSHIDA_H__


/*
 * jth : if N_V, exp(V_jth) 
 *       if N_N, exp(W_jth) where W_j=0.5*s*V_0+sqrt(s)*Z^i_j*V_i
 */
void sde_wa_yoshida(SDE_WA_SLTN *X, double s, int jth);

void exp_2th(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length);
void exp_2th_Kai(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length);
void exp_2th_Omega(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length);
void exp_2th_c3(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length);
void exp_2th_Heston(SDE_WA_SLTN *X, double s, int jth, int d, double w, double W[], int length);

#endif
