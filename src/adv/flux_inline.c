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
/* Functions to be used inline with multi-dimensional flux limiter.  */
                                  /*------------------------------------------*/
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
                                  /*------------------------------------------*/
static inline void get_maxmin_us_2d (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int fid, PREC * u_face, PREC * sc_old, PREC * sc_min,
  PREC * sc_max)
{
int nbrid;
int offset = 1;

if (u_face[fid]*nrm > LT(0.)) offset = 0;

nbrid = ( iY+(offset)*Y_nbr )*nX
        + iX+(offset)*X_nbr;
sc_min[fid] = *(sc_old+nbrid);
sc_max[fid] = *(sc_old+nbrid);

nbrid = ( iY+(offset)*Y_nbr + (1-Y_nbr*nrm) )*nX
        + iX+(offset)*X_nbr + (1-X_nbr*nrm);
sc_min[fid] = MIN(sc_min[fid], *(sc_old+nbrid));
sc_max[fid] = MAX(sc_max[fid], *(sc_old+nbrid));

nbrid = ( iY+(offset)*Y_nbr - (1-Y_nbr*nrm) )*nX
        + iX+(offset)*X_nbr - (1-X_nbr*nrm);
sc_min[fid] = MIN(sc_min[fid], *(sc_old+nbrid));
sc_max[fid] = MAX(sc_max[fid], *(sc_old+nbrid));
}
                                  /*------------------------------------------*/
static inline void limit_inflow_2d (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int nrm, int fid, int tid, PREC * u_face, PREC delta_t, PREC * sc_old,
  PREC * sc_min,
  PREC * sc_max, PREC * sc_face, int * out, PREC * outflow, PREC * inflow,
  PREC * influx_min, PREC * influx_max, PREC * abs_min, PREC * abs_max)
{
int nbrid;
PREC crntxnrm;
crntxnrm = u_face[fid]*delta_t*DELTA_X_INV*nrm;
if (crntxnrm > 0) {
    out[fid] = 1;
    *outflow += crntxnrm;
    nbrid = ( iY+Y_nbr )*nX + iX+X_nbr;
    sc_face[fid] = MIN(sc_face[fid], MAX(*(sc_old+nbrid), sc_max[fid]));
    sc_face[fid] = MAX(sc_face[fid], MIN(*(sc_old+nbrid), sc_min[fid]));
}
else {
    *inflow -= crntxnrm;
    sc_face[fid] = MIN(sc_face[fid], MAX(*(sc_old+tid), sc_max[fid]));
    sc_face[fid] = MAX(sc_face[fid], MIN(*(sc_old+tid), sc_min[fid]));
    sc_min[fid] = MIN(sc_min[fid], sc_face[fid]);
    sc_max[fid] = MAX(sc_max[fid], sc_face[fid]);
    *influx_min -= crntxnrm*sc_min[fid];
    *influx_max -= crntxnrm*sc_max[fid];
    *abs_max = MAX(*abs_max, sc_max[fid]);
    *abs_min = MIN(*abs_min, sc_min[fid]);
}
}
                                  /*------------------------------------------*/
static inline void write_fv (int iX, int iY, int nX, int nY, int X_nbr,
  int Y_nbr, int fv_offset, int fid, PREC * u_face, PREC * sc_fv, PREC * sc_face,
  PREC * sc_face_grad, PREC alpha, PREC out_min, PREC out_max, int * out)
{
if (out[fid]) {
    sc_face[fid] = MIN(sc_face[fid], out_max);
    sc_face[fid] = MAX(sc_face[fid], out_min);
    *(sc_fv + (iY+(Y_nbr))*nX + iX+(X_nbr) + nX*nY*(fv_offset))
      = DELTA_X_INV*alpha*sc_face_grad[fid] - sc_face[fid]*u_face[fid];
}
}
                                  /*------------------------------------------*/

