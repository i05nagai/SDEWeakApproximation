/*
 * sde_wa_butcher.h
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
#include "sde_wa.h"

#ifndef __SDE_WA_BUTCHER_H__
#define __SDE_WA_BUTCHER_H__

/*
 * jth : if N_V, exp(V_jth) 
 *       if N_N, exp(W_jth) where W_j=0.5*s*V_0+sqrt(s)*Z^i_j*V_i
 */
void sde_wa_butcher(SDE_WA_SLTN *X, double s, int jth);

#endif
