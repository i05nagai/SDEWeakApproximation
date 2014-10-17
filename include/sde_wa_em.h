/*
 * sde_wa_em.h 
 */

/*
 * $Author$
 * $Revision$
 * $Source$
 * $Date$
 *
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



#ifndef _SDE_WA_EM_H_
#define _SDE_WA_EM_H_

#include "sde_wa.h"

void Strt_to_Ito_drift(SDE_WA_SLTN *X);
int sde_wa_em(SDE_WA_SLTN *X, double s);

#endif
