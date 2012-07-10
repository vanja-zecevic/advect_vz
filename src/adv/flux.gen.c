/* Copyright (C) 2012 Vanja Zecevic
   Contact vanja.zecevic@sydney.uni.edu.au

   This file is part of advect_vz

   advect_vz is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   advect_vz is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include "src/adv/advection.h"
#include "src/adv/core_face.c"
#include "src/adv/flux_inline.c"
#include "src/adv/flux.h"
#include "src/adv/core.h"
#include "src/adv/bc.h"
/*----------------------------------------------------------------------------*/
/* Generated main loops.  */
                                  /*------------------------------------------*/
#pragma vzg split SCHEME [up laxwend quickest_a quickest_b upmulti utopia \
  utopia_simp]

void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_fvflux_2D)
    (phi_a, phi_fv, u_buff, iX, iY, nX, nY, delta_t, alpha);

if (bc_phi == 1)
    flux_2D_BC_1(phi_fv, u_buff, nX, nY, DELTA_X);
else if (bc_phi == 2)
    facevals_2D_BC_2 (phi_a, phi_fv, nX, nY, DELTA_X);

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(gather_fv_2D)
    (phi_a, phi_fv, iX, iY, nX, nY, delta_t);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_fvflux_2D) (PREC * sc_old, PREC * sc_fv,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;
int out[4]={0, 0, 0, 0};
PREC sc_face[4];
PREC sc_face_grad[4];
PREC sc_min[4];
PREC sc_max[4];
PREC u_face[4];
PREC u_face_t = LT(0.);
PREC inflow = LT(0.);
PREC outflow = LT(0.);
PREC influx_min = LT(0.);
PREC influx_max = LT(0.);
PREC abs_max = -INFINITY;
PREC abs_min = INFINITY;
PREC out_max;
PREC out_min;
PREC flowchange;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, (u_face+0),   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  *(u_face+0), u_face_t, (sc_face+0), sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, (sc_face_grad+0), sc_old);
get_maxmin_us_2d         (iX, iY, nX, nY, 1, 0, 1, 0, u_face,
  sc_old, sc_min, sc_max); 
limit_inflow_2d          (iX, iY, nX, nY, 1, 0, 1, 0, tid, u_face, delta_t,
  sc_old, sc_min, sc_max, sc_face, out, &outflow, &inflow, &influx_min,
  &influx_max, &abs_min, &abs_max);

get_uface_2d             (iX, iY, nX, nY,-1, 0, 0, tid, (u_face+2), u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  *(u_face+2), u_face_t, (sc_face+2), sc_old); 
get_face_grad_2d         (iX, iY, nX, nY,-1, 0, tid, (sc_face_grad+2), sc_old);
get_maxmin_us_2d         (iX, iY, nX, nY,-1, 0,-1, 2, u_face,
  sc_old, sc_min, sc_max); 
limit_inflow_2d          (iX, iY, nX, nY,-1, 0,-1, 2, tid, u_face, delta_t,
  sc_old, sc_min, sc_max, sc_face, out, &outflow, &inflow, &influx_min,
  &influx_max, &abs_min, &abs_max);

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, (u_face+1), u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  *(u_face+1), u_face_t, (sc_face+1), sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, (sc_face_grad+1), sc_old);
get_maxmin_us_2d         (iX, iY, nX, nY, 0, 1, 1, 1, u_face,
  sc_old, sc_min, sc_max); 
limit_inflow_2d          (iX, iY, nX, nY, 0, 1, 1, 1, tid, u_face, delta_t,
  sc_old, sc_min, sc_max, sc_face, out, &outflow, &inflow, &influx_min,
  &influx_max, &abs_min, &abs_max);

get_uface_2d             (iX, iY, nX, nY, 0,-1, 1, tid, (u_face+3), u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  *(u_face+3), u_face_t, (sc_face+3), sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0,-1, tid, (sc_face_grad+3), sc_old);
get_maxmin_us_2d         (iX, iY, nX, nY, 0,-1,-1, 3, u_face,
  sc_old, sc_min, sc_max);
limit_inflow_2d          (iX, iY, nX, nY, 0,-1,-1, 3, tid, u_face, delta_t,
  sc_old, sc_min, sc_max, sc_face, out, &outflow, &inflow, &influx_min,
  &influx_max, &abs_min, &abs_max);

flowchange = LT(1.) + inflow - outflow;
abs_max = MAX(abs_max, *(sc_old+tid));
abs_min = MIN(abs_min, *(sc_old+tid));
out_min = (*(sc_old+tid) + influx_max - abs_max*flowchange)/outflow;
out_max = (*(sc_old+tid) + influx_min - abs_min*flowchange)/outflow;

/* Each cell can write out a out flow face value.  */
write_fv (iX, iY, nX, nY, 0, 0, 0, 0, u_face, sc_fv, sc_face, sc_face_grad,
  alpha, out_min, out_max, out);
write_fv (iX, iY, nX, nY,-1, 0, 0, 2, u_face, sc_fv, sc_face, sc_face_grad,
  alpha, out_min, out_max, out);
write_fv (iX, iY, nX, nY, 0, 0, 1, 1, u_face, sc_fv, sc_face, sc_face_grad,
  alpha, out_min, out_max, out);
write_fv (iX, iY, nX, nY, 0,-1, 1, 3, u_face, sc_fv, sc_face, sc_face_grad,
  alpha, out_min, out_max, out);

}
#pragma vzg splitend

