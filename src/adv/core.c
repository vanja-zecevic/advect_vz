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
#include "src/adv/core.h"
#include "src/adv/bc.h"
/*----------------------------------------------------------------------------*/
/* Generated main loops.  */
                                  /*------------------------------------------*/
#define SCHEME up
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_up
#define SCHEME_NUM 0
void APND_SCHEME(update_2D_nrm) (PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

PREC * phi_old;
PREC * phi_new;

if (iT%2==0) {
    phi_old = phi_a;
    phi_new = phi_b;
}
else {
    phi_old = phi_b;
    phi_new = phi_a;
}

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_cell_2D)
    (phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t, alpha);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_cell_2D) (PREC * sc_old, PREC * sc_new,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC sum_face_a = LT(0.);
PREC sum_face_d = LT(0.);
PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void APND_SCHEME(update_2D_face) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_fv_2D)
    (phi_a, phi_fv, u_buff, iX, iY, nX, nY, delta_t, alpha);

if (bc_phi == 1)
    facevals_2D_BC_1(phi_fv, nX, nY, DELTA_X);
else if (bc_phi == 2)
    facevals_2D_BC_2(phi_a, phi_fv, nX, nY, DELTA_X);

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  gather_fv_2D
    (phi_a, phi_fv, iX, iY, nX, nY, delta_t);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_fv_2D) (PREC * sc_old, PREC * sc_fv,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
*(sc_fv + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
*(sc_fv + nX*nY + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

}

#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME laxwend
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_laxwend
#define SCHEME_NUM 1
void APND_SCHEME(update_2D_nrm) (PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

PREC * phi_old;
PREC * phi_new;

if (iT%2==0) {
    phi_old = phi_a;
    phi_new = phi_b;
}
else {
    phi_old = phi_b;
    phi_new = phi_a;
}

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_cell_2D)
    (phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t, alpha);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_cell_2D) (PREC * sc_old, PREC * sc_new,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC sum_face_a = LT(0.);
PREC sum_face_d = LT(0.);
PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void APND_SCHEME(update_2D_face) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_fv_2D)
    (phi_a, phi_fv, u_buff, iX, iY, nX, nY, delta_t, alpha);

if (bc_phi == 1)
    facevals_2D_BC_1(phi_fv, nX, nY, DELTA_X);
else if (bc_phi == 2)
    facevals_2D_BC_2(phi_a, phi_fv, nX, nY, DELTA_X);

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  gather_fv_2D
    (phi_a, phi_fv, iX, iY, nX, nY, delta_t);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_fv_2D) (PREC * sc_old, PREC * sc_fv,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
*(sc_fv + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
*(sc_fv + nX*nY + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

}

#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME quickest_a
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_quickest_a
#define SCHEME_NUM 2
void APND_SCHEME(update_2D_nrm) (PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

PREC * phi_old;
PREC * phi_new;

if (iT%2==0) {
    phi_old = phi_a;
    phi_new = phi_b;
}
else {
    phi_old = phi_b;
    phi_new = phi_a;
}

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_cell_2D)
    (phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t, alpha);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_cell_2D) (PREC * sc_old, PREC * sc_new,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC sum_face_a = LT(0.);
PREC sum_face_d = LT(0.);
PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void APND_SCHEME(update_2D_face) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_fv_2D)
    (phi_a, phi_fv, u_buff, iX, iY, nX, nY, delta_t, alpha);

if (bc_phi == 1)
    facevals_2D_BC_1(phi_fv, nX, nY, DELTA_X);
else if (bc_phi == 2)
    facevals_2D_BC_2(phi_a, phi_fv, nX, nY, DELTA_X);

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  gather_fv_2D
    (phi_a, phi_fv, iX, iY, nX, nY, delta_t);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_fv_2D) (PREC * sc_old, PREC * sc_fv,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
*(sc_fv + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
*(sc_fv + nX*nY + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

}

#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME quickest_b
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_quickest_b
#define SCHEME_NUM 3
void APND_SCHEME(update_2D_nrm) (PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

PREC * phi_old;
PREC * phi_new;

if (iT%2==0) {
    phi_old = phi_a;
    phi_new = phi_b;
}
else {
    phi_old = phi_b;
    phi_new = phi_a;
}

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_cell_2D)
    (phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t, alpha);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_cell_2D) (PREC * sc_old, PREC * sc_new,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC sum_face_a = LT(0.);
PREC sum_face_d = LT(0.);
PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void APND_SCHEME(update_2D_face) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_fv_2D)
    (phi_a, phi_fv, u_buff, iX, iY, nX, nY, delta_t, alpha);

if (bc_phi == 1)
    facevals_2D_BC_1(phi_fv, nX, nY, DELTA_X);
else if (bc_phi == 2)
    facevals_2D_BC_2(phi_a, phi_fv, nX, nY, DELTA_X);

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  gather_fv_2D
    (phi_a, phi_fv, iX, iY, nX, nY, delta_t);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_fv_2D) (PREC * sc_old, PREC * sc_fv,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
*(sc_fv + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
*(sc_fv + nX*nY + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

}

#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME upmulti
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_upmulti
#define SCHEME_NUM 4
void APND_SCHEME(update_2D_nrm) (PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

PREC * phi_old;
PREC * phi_new;

if (iT%2==0) {
    phi_old = phi_a;
    phi_new = phi_b;
}
else {
    phi_old = phi_b;
    phi_new = phi_a;
}

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_cell_2D)
    (phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t, alpha);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_cell_2D) (PREC * sc_old, PREC * sc_new,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC sum_face_a = LT(0.);
PREC sum_face_d = LT(0.);
PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void APND_SCHEME(update_2D_face) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_fv_2D)
    (phi_a, phi_fv, u_buff, iX, iY, nX, nY, delta_t, alpha);

if (bc_phi == 1)
    facevals_2D_BC_1(phi_fv, nX, nY, DELTA_X);
else if (bc_phi == 2)
    facevals_2D_BC_2(phi_a, phi_fv, nX, nY, DELTA_X);

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  gather_fv_2D
    (phi_a, phi_fv, iX, iY, nX, nY, delta_t);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_fv_2D) (PREC * sc_old, PREC * sc_fv,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
*(sc_fv + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
*(sc_fv + nX*nY + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

}

#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME utopia
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_utopia
#define SCHEME_NUM 5
void APND_SCHEME(update_2D_nrm) (PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

PREC * phi_old;
PREC * phi_new;

if (iT%2==0) {
    phi_old = phi_a;
    phi_new = phi_b;
}
else {
    phi_old = phi_b;
    phi_new = phi_a;
}

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_cell_2D)
    (phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t, alpha);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_cell_2D) (PREC * sc_old, PREC * sc_new,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC sum_face_a = LT(0.);
PREC sum_face_d = LT(0.);
PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void APND_SCHEME(update_2D_face) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_fv_2D)
    (phi_a, phi_fv, u_buff, iX, iY, nX, nY, delta_t, alpha);

if (bc_phi == 1)
    facevals_2D_BC_1(phi_fv, nX, nY, DELTA_X);
else if (bc_phi == 2)
    facevals_2D_BC_2(phi_a, phi_fv, nX, nY, DELTA_X);

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  gather_fv_2D
    (phi_a, phi_fv, iX, iY, nX, nY, delta_t);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_fv_2D) (PREC * sc_old, PREC * sc_fv,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
*(sc_fv + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
*(sc_fv + nX*nY + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

}

#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME utopia_simp
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_utopia_simp
#define SCHEME_NUM 6
void APND_SCHEME(update_2D_nrm) (PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

PREC * phi_old;
PREC * phi_new;

if (iT%2==0) {
    phi_old = phi_a;
    phi_new = phi_b;
}
else {
    phi_old = phi_b;
    phi_new = phi_a;
}

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_cell_2D)
    (phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t, alpha);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_cell_2D) (PREC * sc_old, PREC * sc_new,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC sum_face_a = LT(0.);
PREC sum_face_d = LT(0.);
PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d             (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void APND_SCHEME(update_2D_face) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  APND_SCHEME(update_fv_2D)
    (phi_a, phi_fv, u_buff, iX, iY, nX, nY, delta_t, alpha);

if (bc_phi == 1)
    facevals_2D_BC_1(phi_fv, nX, nY, DELTA_X);
else if (bc_phi == 2)
    facevals_2D_BC_2(phi_a, phi_fv, nX, nY, DELTA_X);

#pragma omp parallel for private(iX, iY)
for (iX=2; iX<(nX-2); iX++) for (iY=2; iY<(nY-2); iY++)
  gather_fv_2D
    (phi_a, phi_fv, iX, iY, nX, nY, delta_t);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void APND_SCHEME(update_fv_2D) (PREC * sc_old, PREC * sc_fv,
  PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha)
{
int tid;

PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

/* x dir */
get_uface_2d             (iX, iY, nX, nY, 1, 0, 0, tid, &u_face,   u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
*(sc_fv + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

/* y dir */
get_uface_2d             (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
#if SCHEME_NUM == 4 || SCHEME_NUM == 5 || SCHEME_NUM == 6
get_uface_2d             (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
#endif
APND_SCHEME(get_face_2d) (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d         (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
*(sc_fv + nX*nY + tid) = DELTA_X_INV*alpha*sc_face_grad - sc_face*u_face;

}

#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM


/*----------------------------------------------------------------------------*/
inline void gather_fv_2D (PREC * sc_old, PREC * sc_fv,
  int iX, int iY, int nX, int nY, PREC delta_t)
{
int tid;

PREC sum_face = LT(0.);
PREC sc_face;

tid = iY*nX + iX;

/* x dir */
read_face (iX, iY, nX, nY, 0, 0, 0, &sc_face, sc_fv);
sum_face += sc_face;

read_face (iX, iY, nX, nY,-1, 0, 0, &sc_face, sc_fv);
sum_face -= sc_face;

/* y dir */
read_face (iX, iY, nX, nY, 0, 0, 1, &sc_face, sc_fv);
sum_face += sc_face;

read_face (iX, iY, nX, nY, 0,-1, 1, &sc_face, sc_fv);
sum_face -= sc_face;

*(sc_old + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t*sum_face;

}

