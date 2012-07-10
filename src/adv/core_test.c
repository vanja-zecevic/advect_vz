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
#include "src/adv/core_test.h"
#include "src/adv/bc.h"

/*----------------------------------------------------------------------------*/
/* Lliffe main loop.  */
                                  /*------------------------------------------*/
void update_2D_quickest_a_lliffe (PREC ** phi_a, PREC ** phi_b, 
  PREC * u_buff, int iT, int nX,
  int nY, int bc_phi, PREC delta_t, PREC alpha)
{
int iX;
int iY;

PREC ** phi_old;
PREC ** phi_new;

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
  update_cell_2D_quickest_a_lliffe
    (phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t, alpha);

if      (bc_phi == 1)
    Phi_2D_BC_1(*phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(*phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
inline void update_cell_2D_quickest_a_lliffe (PREC ** sc_old, PREC ** sc_new,
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
get_uface_2d                 (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
get_face_2d_quickest_a_lliffe(iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d_lliffe      (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d                 (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
get_face_2d_quickest_a_lliffe(iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d_lliffe      (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d                 (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
get_face_2d_quickest_a_lliffe(iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d_lliffe      (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d                 (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
get_face_2d_quickest_a_lliffe(iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d_lliffe      (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

sc_new[iY][iX] = sc_old[iY][iX] + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
/* Branching main loops.  */
                                  /*------------------------------------------*/
/*----------------------------------------------------------------------------*/
void update_2D_BRANCH(PREC * phi_a, PREC * phi_b, PREC * u_buff, int iT,
  int nX, int nY, int bc_phi, PREC delta_t, PREC alpha, int scheme)
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
  update_cell_2D_BRANCH(phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t,
  alpha, scheme);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
void update_cell_2D_BRANCH(PREC * sc_old, PREC * sc_new, PREC * u_buff,
  int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha, int scheme)
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
get_uface_2d            (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
if (scheme==0)
  get_face_2d_up        (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==1)
  get_face_2d_laxwend   (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==2)
  get_face_2d_quickest_a(iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==3)
  get_face_2d_quickest_b(iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==4) {
    get_uface_2d           (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_upmulti    (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else if (scheme==5) {
    get_uface_2d           (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_utopia_simp(iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else if (scheme==6) {
    get_uface_2d           (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_utopia     (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
} 
else sc_face = LT(0.);
get_face_grad_2d        (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d            (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
if (scheme==0)
  get_face_2d_up        (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==1)
  get_face_2d_laxwend   (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==2)
  get_face_2d_quickest_a(iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==3)
  get_face_2d_quickest_b(iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==4) {
    get_uface_2d           (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_upmulti    (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else if (scheme==5) {
    get_uface_2d           (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_utopia_simp(iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else if (scheme==6) {
    get_uface_2d           (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_utopia     (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else sc_face = LT(0.);
get_face_grad_2d        (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d            (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
if (scheme==0)
  get_face_2d_up        (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==1)
  get_face_2d_laxwend   (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==2)
  get_face_2d_quickest_a(iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==3)
  get_face_2d_quickest_b(iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==4) {
    get_uface_2d           (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
    get_face_2d_upmulti    (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else if (scheme==5) {
    get_uface_2d           (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
    get_face_2d_utopia_simp(iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else if (scheme==6) {
    get_uface_2d           (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
    get_face_2d_utopia     (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else sc_face = LT(0.);
get_face_grad_2d        (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d            (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
if (scheme==0)
  get_face_2d_up        (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==1)
  get_face_2d_laxwend   (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==2)
  get_face_2d_quickest_a(iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==3)
  get_face_2d_quickest_b(iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
    u_face, u_face_t, &sc_face, sc_old); 
else if (scheme==4) {
    get_uface_2d           (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
    get_face_2d_upmulti    (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else if (scheme==5) {
    get_uface_2d           (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
    get_face_2d_utopia_simp(iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else if (scheme==6) {
    get_uface_2d           (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
    get_face_2d_utopia     (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old);
}
else sc_face = LT(0.);
get_face_grad_2d        (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void update_2D_BRANCH2(PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha, int scheme)
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
  update_cell_2D_BRANCH2(phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t,
  alpha, scheme);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
void update_cell_2D_BRANCH2(PREC * sc_old, PREC * sc_new, PREC * u_buff,
  int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha, int scheme)
{
int tid;

PREC sum_face_a = LT(0.);
PREC sum_face_d = LT(0.);
PREC u_face;
PREC u_face_t = LT(0.);
PREC sc_face;
PREC sc_face_grad;

tid = iY*nX + iX;

if (scheme==0) {
    /* x dir */
    get_uface_2d    (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
    get_face_2d_up  (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d(iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d    (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
    get_face_2d_up  (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d(iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;

    /* y dir */
    get_uface_2d    (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
    get_face_2d_up  (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d(iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d    (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
    get_face_2d_up  (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d(iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;
}
else if (scheme==1) {
    /* x dir */
    get_uface_2d       (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
    get_face_2d_laxwend(iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d   (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d       (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
    get_face_2d_laxwend(iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d   (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;

    /* y dir */
    get_uface_2d       (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
    get_face_2d_laxwend(iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d   (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d       (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
    get_face_2d_laxwend(iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d   (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;
}
else if (scheme==2) {
    /* x dir */
    get_uface_2d          (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
    get_face_2d_quickest_a(iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d          (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
    get_face_2d_quickest_a(iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;

    /* y dir */
    get_uface_2d          (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
    get_face_2d_quickest_a(iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d          (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
    get_face_2d_quickest_a(iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;
}
else if (scheme==3) {
    /* x dir */
    get_uface_2d          (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
    get_face_2d_quickest_b(iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d          (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
    get_face_2d_quickest_b(iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;

    /* y dir */
    get_uface_2d          (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
    get_face_2d_quickest_b(iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d          (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
    get_face_2d_quickest_b(iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;
}
else if (scheme==4) {
    /* x dir */
    get_uface_2d          (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
    get_uface_2d          (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_upmulti   (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d          (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
    get_uface_2d          (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_upmulti   (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;

    /* y dir */
    get_uface_2d          (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
    get_uface_2d          (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
    get_face_2d_upmulti   (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d          (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
    get_uface_2d          (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
    get_face_2d_upmulti   (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;
}
else if (scheme==5) {
    /* x dir */
    get_uface_2d           (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
    get_uface_2d           (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_utopia_simp(iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d           (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
    get_uface_2d           (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_utopia_simp(iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d       (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;

    /* y dir */
    get_uface_2d           (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
    get_uface_2d           (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
    get_face_2d_utopia_simp(iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d       (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d           (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
    get_uface_2d           (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
    get_face_2d_utopia_simp(iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d       (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;
}
else if (scheme==6) {
    /* x dir */
    get_uface_2d           (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
    get_uface_2d           (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_utopia     (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d      (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d           (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
    get_uface_2d           (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
    get_face_2d_utopia     (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d       (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;

    /* y dir */
    get_uface_2d           (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
    get_uface_2d           (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
    get_face_2d_utopia     (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d       (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
    sum_face_a += sc_face*u_face;
    sum_face_d += sc_face_grad;

    get_uface_2d           (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
    get_uface_2d           (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
    get_face_2d_utopia     (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
      u_face, u_face_t, &sc_face, sc_old); 
    get_face_grad_2d       (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
    sum_face_a -= sc_face*u_face;
    sum_face_d += sc_face_grad;
}
else sc_face = LT(0.);

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/
void update_2D_PTRS(PREC * phi_a, PREC * phi_b, PREC * u_buff,
  int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha, int scheme)
{
int iX;
int iY;

PREC * phi_old;
PREC * phi_new;

void  (*get_face_2d_ptr) (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC dummy,
  PREC * sc_face, PREC * sc_old);

if      (scheme==0)
  get_face_2d_ptr = &get_face_2d_up; 
else if (scheme==1)
  get_face_2d_ptr = &get_face_2d_laxwend; 
else if (scheme==2)
  get_face_2d_ptr = &get_face_2d_quickest_a; 
else if (scheme==3)
  get_face_2d_ptr = &get_face_2d_quickest_b; 
else if (scheme==4)
  get_face_2d_ptr = &get_face_2d_upmulti; 
else if (scheme==5)
  get_face_2d_ptr = &get_face_2d_utopia_simp; 
else if (scheme==6)
  get_face_2d_ptr = &get_face_2d_utopia; 
else exit(1);

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
  update_cell_2D_PTRS(phi_old, phi_new, u_buff, iX, iY, nX, nY, delta_t,
  alpha, scheme, get_face_2d_ptr);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_new, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_new, nX, nY, DELTA_X);

}
/*----------------------------------------------------------------------------*/
void update_cell_2D_PTRS(PREC * sc_old, PREC * sc_new, PREC * u_buff,
  int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha, int scheme,
  void (*get_face_2d_ptr) (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC dummy,
  PREC * sc_face, PREC * sc_old))
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
get_uface_2d     (iX, iY, nX, nY, 1, 0, 0, tid, &u_face, u_buff);
if (scheme==4 || scheme==5 || scheme==6)
  get_uface_2d   (iX, iY, nX, nY, 1, 0, 1, tid, &u_face_t, u_buff);
get_face_2d_ptr  (iX, iY, nX, nY, 1, 0, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d (iX, iY, nX, nY, 1, 0, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d     (iX, iY, nX, nY,-1, 0, 0, tid, &u_face, u_buff);
if (scheme==4 || scheme==5 || scheme==6)
  get_uface_2d   (iX, iY, nX, nY,-1, 0, 1, tid, &u_face_t, u_buff);
get_face_2d_ptr  (iX, iY, nX, nY,-1, 0,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d (iX, iY, nX, nY,-1, 0, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

/* y dir */
get_uface_2d     (iX, iY, nX, nY, 0, 1, 1, tid, &u_face, u_buff);
if (scheme==4 || scheme==5 || scheme==6)
  get_uface_2d   (iX, iY, nX, nY, 0, 1, 0, tid, &u_face_t, u_buff);
get_face_2d_ptr  (iX, iY, nX, nY, 0, 1, 1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d (iX, iY, nX, nY, 0, 1, tid, &sc_face_grad, sc_old);
sum_face_a += sc_face*u_face;
sum_face_d += sc_face_grad;

get_uface_2d     (iX, iY, nX, nY, 0,-1, 1, tid, &u_face, u_buff);
if (scheme==4 || scheme==5 || scheme==6)
  get_uface_2d   (iX, iY, nX, nY, 0,-1, 0, tid, &u_face_t, u_buff);
get_face_2d_ptr  (iX, iY, nX, nY, 0,-1,-1, tid, delta_t,
  u_face, u_face_t, &sc_face, sc_old); 
get_face_grad_2d (iX, iY, nX, nY, 0,-1, tid, &sc_face_grad, sc_old);
sum_face_a -= sc_face*u_face;
sum_face_d += sc_face_grad;

*(sc_new + tid) = *(sc_old + tid) + DELTA_X_INV*delta_t
      *(DELTA_X_INV*alpha*sum_face_d - sum_face_a);

}
/*----------------------------------------------------------------------------*/

