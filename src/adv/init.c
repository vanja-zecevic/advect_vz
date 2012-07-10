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
#include "src/adv/init.h"

void init_u_rigid(PREC * u_buff, int nX, int nY, int nZ, PREC u_init_max)
{
/* Initializes the velocity field to rigid body rotation about an axis paralell
   to the z-axis in the center of the domain. This works in 2D as well if
   nZ = 1.  */
int iX;
int iY;
int iZ;

PREC x_zero;
PREC y_zero;
PREC r_x;
PREC r_y;
PREC omega;

omega = u_init_max/( (PREC)(nX-4)*LT(0.5) );
x_zero = (PREC)(nX-4)*LT(0.5);
y_zero = (PREC)(nY-4)*LT(0.5);

for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++) {
    r_x = iX - LT(1.5) - x_zero;
    r_y = iY - LT(1.5) - y_zero;
    *(u_buff + 0*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) = -r_y*omega;
    *(u_buff + 1*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) =  r_x*omega;
    *(u_buff + 2*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) =  LT(0.);
}
}
/*----------------------------------------------------------------------------*/
void init_u_const(PREC * u_buff, int nX, int nY, int nZ, PREC u_init_max,
  PREC u_init_angle)
{
/* Constant velocity initialization.  */
int iX;
int iY;
int iZ;

for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++) {
    *(u_buff + 0*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) =
      u_init_max*cos(u_init_angle);
    *(u_buff + 1*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) =
      u_init_max*sin(u_init_angle);
    *(u_buff + 2*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) = LT(0.);
}
}
/*----------------------------------------------------------------------------*/
void init_u_tvort(PREC * u_buff, int nX, int nY, int nZ, PREC u_init_max)
{
int iX;
int iY;
int iZ;
PREC kx;
PREC ky;

kx = LT(2.0)*VZ_PI/(PREC)(nX-4);
ky = LT(2.0)*VZ_PI/(PREC)(nY-4);

for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++) {
    *(u_buff + 0*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) =
       u_init_max*sin( kx*(iX-LT(1.5)) )*cos( ky*(iY-LT(1.5)) );
    *(u_buff + 1*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) =
      -u_init_max*cos( kx*(iX-LT(1.5)) )*sin( ky*(iY-LT(1.5)) );
    *(u_buff + 2*nX*nY*nZ + iZ*nX*nY + iY*nX + iX) =  LT(0.);
}
}
/*----------------------------------------------------------------------------*/
void init_phi_gauss(PREC * phi, int nX, int nY, int nZ, PREC x_zero,
  PREC y_zero, PREC z_zero, PREC phi_init_max, PREC B)
{
int iX;
int iY;
int iZ;

PREC delta_x;
PREC delta_y;
PREC delta_z;

PREC twoBsq;

twoBsq = -2*B*B;

for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++) {
    delta_x = iX - LT(1.5) - x_zero;
    delta_y = iY - LT(1.5) - y_zero;
    delta_z = iZ - LT(1.5) - z_zero;
    *(phi + iZ*nX*nY + iY*nX + iX) = phi_init_max*
      exp((delta_x*delta_x + delta_y*delta_y + delta_z*delta_z)/twoBsq);
}

}
/*----------------------------------------------------------------------------*/
void init_phi_cyl(PREC * phi, int nX, int nY, int nZ, PREC x_zero,
  PREC y_zero, PREC z_zero, PREC phi_init_max, PREC radius)
{
int iX;
int iY;
int iZ;

PREC delta_x;
PREC delta_y;

for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++) {
    delta_x = iX - LT(1.5) - x_zero;
    delta_y = iY - LT(1.5) - y_zero;
    if ( sqrt(delta_x*delta_x + delta_y*delta_y) <= radius)
      *(phi + iZ*nX*nY + iY*nX + iX) = phi_init_max;
    else
      *(phi + iZ*nX*nY + iY*nX + iX) = LT(0.);

}

}
/*----------------------------------------------------------------------------*/
void init_phi_step(PREC * phi, int nX, int nY, int nZ, PREC phi_init_max,
  PREC h_step)
{
int iX;
int iY;
int iZ;

for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++) {
    *(phi + iZ*nX*nY + iY*nX + iX) = phi_init_max*( LT(2.) +
      tanh( (LT(0.25) - (PREC)(iX-LT(1.5))/(PREC)(nX-4))/h_step ) - 
      tanh( (LT(0.75) - (PREC)(iX-LT(1.5))/(PREC)(nX-4))/h_step ));
}

}
/*----------------------------------------------------------------------------*/
void init_phi_const(PREC * phi, int nX, int nY, int nZ, PREC phi_init_max)
{
int iX;
int iY;
int iZ;

for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++) {
    *(phi + iZ*nX*nY + iY*nX + iX) = phi_init_max;
}
}
/*----------------------------------------------------------------------------*/

