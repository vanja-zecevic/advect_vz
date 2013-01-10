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
/*----------------------------------------------------------------------------*/
/*  Functions to be used inline.  */
                                  /*------------------------------------------*/
static inline void get_face_2d_up (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC dummy,
  PREC * sc_face, PREC * sc_old)
{
int offset = 1;
int nbrid;
PREC windxnrm;

windxnrm = u_face*nrm;
if (windxnrm > LT(0.)) offset = 0;

nbrid = (iY+offset*Y_nbr)*nX + iX+offset*X_nbr;
*sc_face = *(sc_old+nbrid);
}
                                  /*------------------------------------------*/
static inline void get_face_2d_laxwend (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC dummy,
  PREC * sc_face, PREC * sc_old)
{
int nbrid;
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
*sc_face = LT(0.5)*( *(sc_old+nbrid) + *(sc_old+tid)
  - u_face*delta_t*DELTA_X_INV*nrm*( *(sc_old+nbrid) - *(sc_old+tid) ) );
}
                                  /*------------------------------------------*/
static inline void get_face_2d_quickest_a (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC dummy,
  PREC * sc_face, PREC * sc_old)
{
/* quickest_a uses direct memory access however contains more
 * arithmatic (similar speed).  */
int nbrid;
PREC wind = LT(1.);
PREC crnt;
PREC tmp;
PREC sc_face_tmp;

crnt = u_face*delta_t*DELTA_X_INV;
if (crnt < LT(0.)) wind = LT(-1.);
tmp = (LT(1.)-crnt*crnt)*LT(0.08333333333333333);

sc_face_tmp = *(sc_old+tid)*( LT(0.5) + LT(0.5)*crnt*nrm +
  tmp*(LT(1.)+LT(3.)*wind*(nrm)) );
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
sc_face_tmp += *(sc_old+nbrid)*( LT(0.5) - LT(0.5)*crnt*(nrm) +
  tmp*(LT(1.)-LT(3.)*wind*(nrm)) );
nbrid = ( iY-(Y_nbr) )*nX + iX-(X_nbr);
sc_face_tmp += *(sc_old+nbrid)*(
  tmp*(-LT(1.)-(PREC)(wind*(nrm))) );
nbrid = ( iY+2*(Y_nbr) )*nX + iX+2*(X_nbr);
sc_face_tmp += *(sc_old+nbrid)*(
  tmp*(-LT(1.)+(PREC)(wind*(nrm))) );
*sc_face = sc_face_tmp;
}
                                  /*------------------------------------------*/
static inline void get_face_2d_quickest_a_lliffe (int iX, int iY, int nX, int nY,
  int X_nbr, int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC dummy,
  PREC * sc_face, PREC ** sc_old)
{
/* quickest_a uses direct memory access however contains more
 * arithmatic (similar speed).  */
PREC wind = LT(1.);
PREC crnt;
PREC tmp;
PREC sc_face_tmp;

crnt = u_face*delta_t*DELTA_X_INV;

if (crnt < LT(0.)) wind = LT(-1.);

tmp = (LT(1.)-crnt*crnt)*LT(0.08333333333333333);
sc_face_tmp = sc_old[iY][iX]
  *( LT(0.5) + LT(0.5)*crnt*nrm + tmp*(LT(1.)+LT(3.)*wind*(nrm)) );
sc_face_tmp += sc_old[iY+Y_nbr][iX+X_nbr]
  *( LT(0.5) - LT(0.5)*crnt*(nrm) + tmp*(LT(1.)-LT(3.)*wind*(nrm)) );
sc_face_tmp += sc_old[iY-Y_nbr][iX-X_nbr]
  *(tmp*(-LT(1.)-(PREC)(wind*(nrm))));
sc_face_tmp += sc_old[iY+2*Y_nbr][iX+2*X_nbr]
  *(tmp*(-LT(1.)+(PREC)(wind*(nrm))));
*sc_face = sc_face_tmp;
}
                                  /*------------------------------------------*/
static inline void get_face_2d_quickest_b (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC dummy,
  PREC * sc_face, PREC * sc_old)
{
/* quickest_b uses indirect memory access however contains less
 * arithmatic (similar speed).  */
int nbrid;
int offset = 1;
PREC crntxnrm;
PREC tmp;
PREC sc_face_tmp;

crntxnrm = u_face*delta_t*DELTA_X_INV*nrm;
if (crntxnrm > LT(0.)) offset = 0;

/* Lax-Wendroff component.  */
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
sc_face_tmp = LT(0.5)*( *(sc_old+nbrid) + *(sc_old+tid)
  - crntxnrm*( *(sc_old+nbrid) - *(sc_old+tid) ) );

/* Curvature component, X_nbr*nrm is 0 or 1.  */
tmp = (LT(1.)-crntxnrm*crntxnrm)*LT(0.1666666666666666);
nbrid = ( iY+(offset)*Y_nbr )*nX + iX+(offset)*X_nbr;
sc_face_tmp += *(sc_old+nbrid)*LT(2.)*tmp;
nbrid = ( iY+(offset+nrm)*Y_nbr )*nX + iX+(offset+nrm)*X_nbr;
sc_face_tmp -= *(sc_old+nbrid)*tmp;
nbrid = ( iY+(offset-nrm)*Y_nbr )*nX + iX+(offset-nrm)*X_nbr;
sc_face_tmp -= *(sc_old+nbrid)*tmp;

*sc_face = sc_face_tmp;
}
                                  /*------------------------------------------*/
static inline void get_face_2d_utopia_simp (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC u_face_t,
  PREC * sc_face, PREC * sc_old)
{
int nbrid;
int offset = 1;
PREC crntxnrm;
PREC tmp;
PREC sc_face_tmp;

crntxnrm = u_face*delta_t*DELTA_X_INV*nrm;
if (crntxnrm > LT(0.)) offset = 0;

/* Lax-Wendroff component.  */
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
sc_face_tmp = LT(0.5)*( *(sc_old+nbrid) + *(sc_old+tid)
  - crntxnrm*( *(sc_old+nbrid) - *(sc_old+tid) ) );

/* Curvature component, X_nbr*nrm is 0 or 1.  */
tmp = (LT(1.)-crntxnrm*crntxnrm)*LT(0.1666666666666666);
nbrid = ( iY+(offset)*Y_nbr )*nX + iX+(offset)*X_nbr;
sc_face_tmp += *(sc_old+nbrid)*LT(2.)*tmp;
nbrid = ( iY+(offset+nrm)*Y_nbr )*nX + iX+(offset+nrm)*X_nbr;
sc_face_tmp -= *(sc_old+nbrid)*tmp;
nbrid = ( iY+(offset-nrm)*Y_nbr )*nX + iX+(offset-nrm)*X_nbr;
sc_face_tmp -= *(sc_old+nbrid)*tmp;

nbrid = ( iY+(offset)*Y_nbr )*nX
        + iX+(offset)*X_nbr;
sc_face_tmp -= *(sc_old+nbrid)*u_face_t*delta_t*DELTA_X_INV*
  u_face_t*delta_t*DELTA_X_INV/LT(3.);
nbrid = ( iY+(offset)*Y_nbr + (1-Y_nbr*nrm) )*nX
        + iX+(offset)*X_nbr + (1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*u_face_t*delta_t*DELTA_X_INV*(
  - LT(0.25) + u_face_t*delta_t*DELTA_X_INV*LT(0.1666666666666666) );
nbrid = ( iY+(offset)*Y_nbr - (1-Y_nbr*nrm) )*nX
        + iX+(offset)*X_nbr - (1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*u_face_t*delta_t*DELTA_X_INV*(
  LT(0.25) + u_face_t*delta_t*DELTA_X_INV*LT(0.1666666666666666) );

*sc_face = sc_face_tmp;

}
                                  /*------------------------------------------*/
static inline void get_face_2d_utopia (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC u_face_t,
  PREC * sc_face, PREC * sc_old)
{
int nbrid;
int offset = 1;
int iwind_t = 1;

PREC tmp;
PREC tmp_2;
PREC sc_face_tmp;
PREC wind   = LT(1.);
PREC wind_t  = LT(1.);
PREC crnt;
PREC crnt_t;

crnt = u_face*delta_t*DELTA_X_INV;
if (crnt*nrm > LT(0.)) offset = 0;

/* Lax-Wendroff component.  */
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
sc_face_tmp = LT(0.5)*( *(sc_old+nbrid) + *(sc_old+tid)
  - crnt*nrm*( *(sc_old+nbrid) - *(sc_old+tid) ) );

/* Curvature component, X_nbr*nrm is 0 or 1.  */
tmp = (LT(1.)-crnt*crnt)*LT(0.1666666666666666);
nbrid = ( iY+(offset)*Y_nbr )*nX + iX+(offset)*X_nbr;
sc_face_tmp += *(sc_old+nbrid)*LT(2.)*tmp;
nbrid = ( iY+(offset+nrm)*Y_nbr )*nX + iX+(offset+nrm)*X_nbr;
sc_face_tmp -= *(sc_old+nbrid)*tmp;
nbrid = ( iY+(offset-nrm)*Y_nbr )*nX + iX+(offset-nrm)*X_nbr;
sc_face_tmp -= *(sc_old+nbrid)*tmp;

nbrid = ( iY+(offset)*Y_nbr )*nX
        + iX+(offset)*X_nbr;
sc_face_tmp -= *(sc_old+nbrid)*LT(0.5)*u_face_t*delta_t*DELTA_X_INV*
  u_face_t*delta_t*DELTA_X_INV;
nbrid = ( iY+(offset)*Y_nbr + (1-Y_nbr*nrm) )*nX
        + iX+(offset)*X_nbr + (1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*u_face_t*delta_t*DELTA_X_INV*(
  - LT(0.25) + LT(0.25)*u_face_t*delta_t*DELTA_X_INV );
nbrid = ( iY+(offset)*Y_nbr - (1-Y_nbr*nrm) )*nX
        + iX+(offset)*X_nbr - (1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*u_face_t*delta_t*DELTA_X_INV*(
    LT(0.25) + LT(0.25)*u_face_t*delta_t*DELTA_X_INV );

if (crnt   < LT(0.)) wind = LT(-1.);
crnt_t = u_face_t*delta_t*DELTA_X_INV;
if (crnt_t < LT(0.)) {
    wind_t  = LT(-1.);
    iwind_t =    -1;
}
tmp_2 = LT(0.25)*crnt_t*(wind - crnt);
sc_face_tmp += *(sc_old+tid)*wind_t*nrm*tmp_2;
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
sc_face_tmp -= *(sc_old+nbrid)*wind_t*nrm*tmp_2;
nbrid = ( iY-iwind_t*(1-Y_nbr*nrm) )*nX
        + iX-iwind_t*(1-X_nbr*nrm);
sc_face_tmp -= *(sc_old+nbrid)*wind_t*nrm*tmp_2;
nbrid = ( iY+(Y_nbr)-iwind_t*(1-Y_nbr*nrm) )*nX
        + iX+(X_nbr)-iwind_t*(1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*wind_t*nrm*tmp_2;

*sc_face = sc_face_tmp;
}
                                  /*------------------------------------------*/
static inline void get_face_2d_utopia_b (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC u_face_t,
  PREC * sc_face, PREC * sc_old)
{
PREC sc_face_tmp;
PREC wind    = LT(1.);
PREC wind_t  = LT(1.);
PREC crnt;
PREC crnt_t;
PREC tmp_1;
PREC tmp_2;
int iwind   = 1;
int iwind_t = 1;
int nbrid;

crnt   = u_face  *delta_t*DELTA_X_INV;
crnt_t = u_face_t*delta_t*DELTA_X_INV;

if (crnt   < LT(0.)) {
    wind    = LT(-1.);
    iwind   =    -1;
}
if (crnt_t < LT(0.)) {
    wind_t  = LT(-1.);
    iwind_t =    -1;
}
tmp_1 = (LT(1.)-crnt*crnt)*LT(0.0833333333333333333);
tmp_2 = LT(0.25)*crnt_t*(wind - crnt);

/* quickest scheme.  */
sc_face_tmp = *(sc_old+tid)*( LT(0.5) + LT(0.5)*crnt*nrm +
  tmp_1*(LT(1.)+LT(3.)*wind*(nrm)) );
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
sc_face_tmp += *(sc_old+nbrid)*( LT(0.5) - LT(0.5)*crnt*(nrm) +
  tmp_1*(LT(1.)-LT(3.)*wind*(nrm)) );
nbrid = ( iY-(Y_nbr) )*nX + iX-(X_nbr);
sc_face_tmp += *(sc_old+nbrid)*(
  tmp_1*(-LT(1.)-(PREC)(wind*(nrm))) );
nbrid = ( iY+2*(Y_nbr) )*nX + iX+2*(X_nbr);
sc_face_tmp += *(sc_old+nbrid)*(
  tmp_1*(-LT(1.)+(PREC)(wind*(nrm))) );

nbrid = ( iY+(Y_nbr)-Y_nbr*(iwind*nrm+1)/2 )*nX
        + iX+(X_nbr)-X_nbr*(iwind*nrm+1)/2;
sc_face_tmp -= *(sc_old+nbrid)*LT(0.5)*u_face_t*delta_t*DELTA_X_INV*
  u_face_t*delta_t*DELTA_X_INV;
nbrid = ( iY+(Y_nbr)-Y_nbr*(iwind*nrm+1)/2 + (1-Y_nbr*nrm) )*nX
        + iX+(X_nbr)-X_nbr*(iwind*nrm+1)/2 + (1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*u_face_t*delta_t*DELTA_X_INV*(
  - LT(0.25) + LT(0.25)*u_face_t*delta_t*DELTA_X_INV );
nbrid = ( iY+(Y_nbr)-Y_nbr*(iwind*nrm+1)/2 - (1-Y_nbr*nrm) )*nX
        + iX+(X_nbr)-X_nbr*(iwind*nrm+1)/2 - (1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*u_face_t*delta_t*DELTA_X_INV*(
    LT(0.25) + LT(0.25)*u_face_t*delta_t*DELTA_X_INV );

sc_face_tmp += *(sc_old+tid)*wind_t*nrm*tmp_2;
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
sc_face_tmp -= *(sc_old+nbrid)*wind_t*nrm*tmp_2;
nbrid = ( iY-iwind_t*(1-Y_nbr*nrm) )*nX
        + iX-iwind_t*(1-X_nbr*nrm);
sc_face_tmp -= *(sc_old+nbrid)*wind_t*nrm*tmp_2;
nbrid = ( iY+(Y_nbr)-iwind_t*(1-Y_nbr*nrm) )*nX
        + iX+(X_nbr)-iwind_t*(1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*wind_t*nrm*tmp_2;
*sc_face = sc_face_tmp;
}
                                  /*------------------------------------------*/
static inline void get_face_2d_upmulti (int iX, int iY, int nX, int nY,
  int X_nbr, int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face,
  PREC u_face_t, PREC * sc_face, PREC * sc_old)
{
int offset = 1;
int iwind_t = 1;
int nbrid;
PREC wind_t  = LT(1.);
PREC abs_crnt_t;
PREC sc_face_tmp;

if (u_face*nrm > LT(0.)) offset = 0;
if (u_face_t < LT(0.)) {
    wind_t  = LT(-1.);
    iwind_t =    -1;
}
abs_crnt_t = wind_t*u_face_t*delta_t*DELTA_X_INV;

nbrid = ( iY+offset*Y_nbr)*nX
        + iX+offset*X_nbr;
sc_face_tmp = *(sc_old+nbrid)*( LT(1.) - LT(0.5)*abs_crnt_t );
nbrid = ( iY+offset*Y_nbr - iwind_t*(1-Y_nbr*nrm) )*nX
        + iX+offset*X_nbr - iwind_t*(1-X_nbr*nrm);
sc_face_tmp += *(sc_old+nbrid)*( LT(0.5)*abs_crnt_t );

*sc_face = sc_face_tmp;

}
                                  /*------------------------------------------*/
static inline void get_uface_2d (int iX, int iY, int nX, int nY, int X_nbr, int Y_nrb,
  int u_dir, int tid, PREC * uface, PREC * u_buff)
{
int nbrid;
nbrid = ( iY+(Y_nrb) )*nX + iX+(X_nbr);
*uface = LT(0.5)*( *(u_buff+nbrid+u_dir*nX*nY) + *(u_buff+tid+u_dir*nX*nY) );
}
                                  /*------------------------------------------*/
static inline void get_face_grad_2d (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int tid, PREC * sc_face_grad, PREC * sc_old)
{
int nbrid;
nbrid = ( iY+(Y_nbr) )*nX + iX+(X_nbr);
*sc_face_grad = *(sc_old+nbrid) - *(sc_old+tid);
}
                                  /*------------------------------------------*/
static inline void get_face_grad_2d_lliffe (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int tid, PREC * sc_face_grad, PREC ** sc_old)
{
*sc_face_grad = sc_old[iY+Y_nbr][iX+X_nbr] - sc_old[iY][iX];
}
                                  /*------------------------------------------*/
static inline void read_face(int iX, int iY, int nX, int nY, int X_nbr, int Y_nbr,
  int fid, PREC * sc_face, PREC * sc_fv) 
{
*sc_face = *(sc_fv + (iY+(Y_nbr))*nX + iX+(X_nbr) + nX*nY*(fid));
}
                                  /*------------------------------------------*/
static inline void gather_fv_2D (PREC * sc_old, PREC * sc_fv,
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
                                  /*------------------------------------------*/

