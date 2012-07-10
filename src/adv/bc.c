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
#include "src/adv/bc.h"
/*---------------------------------------------------------------------------*/
void Neumann_XZ(PREC * phi, int nX, int nY, int nZ, int Y_loc, int normal,
  PREC grad, PREC delta_x)
{
int iX;
int iZ;
for (iX=0; iX<nX; iX++) for (iZ=0; iZ<nZ; iZ++)
  *(phi + iZ*nX*nY + Y_loc*nX + iX) =
    *(phi + iZ*nX*nY + (Y_loc+normal)*nX + iX)
    - (PREC)normal*delta_x*grad;

}
/*---------------------------------------------------------------------------*/
void Neumann_XY(PREC * phi, int nX, int nY, int nZ, int Z_loc, int normal,
  PREC grad, PREC delta_x)
{
int iX;
int iY;
for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++)
  *(phi + Z_loc*nX*nY + iY*nX + iX) =
    *(phi + (Z_loc+normal)*nX*nY + iY*nX + iX)
    - (PREC)normal*delta_x*grad;

}
/*---------------------------------------------------------------------------*/
void Neumann_YZ(PREC * phi, int nX, int nY, int nZ, int X_loc, int normal,
  PREC grad, PREC delta_x)
{
int iY;
int iZ;
for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++)
  *(phi + iZ*nX*nY + iY*nX + X_loc) =
    *(phi + iZ*nX*nY + iY*nX + (X_loc+normal))
    - (PREC)normal*delta_x*grad;

}
/*---------------------------------------------------------------------------*/
void Zero_XZ(PREC * phi, int nX, int nY, int nZ, int Y_loc)
{
int iX;
int iZ;
for (iX=0; iX<nX; iX++) for (iZ=0; iZ<nZ; iZ++)
  *(phi + iZ*nX*nY + Y_loc*nX + iX) = LT(0.);

}
/*---------------------------------------------------------------------------*/
void Zero_XY(PREC * phi, int nX, int nY, int nZ, int Z_loc)
{
int iX;
int iY;
for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++)
  *(phi + Z_loc*nX*nY + iY*nX + iX) = LT(0.);

}
/*---------------------------------------------------------------------------*/
void Zero_YZ(PREC * phi, int nX, int nY, int nZ, int X_loc)
{
int iY;
int iZ;
for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++)
  *(phi + iZ*nX*nY + iY*nX + X_loc) = LT(0.);

}
/*---------------------------------------------------------------------------*/
void Periodic_XZ_flux(PREC * phi_fv, PREC * u_buff, int nX, int nY, int nZ,
  int Y_low, int Y_high)
{
int iX;
int iZ;
PREC u_face;

for (iX=0; iX<nX; iX++) for (iZ=0; iZ<nZ; iZ++) {
    u_face = LT(0.5)*( *(u_buff + iX + Y_low*nX + iZ*nX*nY + (1)*nX*nY*nZ) 
      + *(u_buff + iX + (Y_low+1)*nX + iZ*nX*nY + (1)*nX*nY*nZ) );
    if (u_face > 0)
      *(phi_fv + iZ*nX*nY + Y_low*nX + iX + (1)*nX*nY*nZ) =
        *(phi_fv + iZ*nX*nY + Y_high*nX + iX + (1)*nX*nY*nZ);
    else
      *(phi_fv + iZ*nX*nY + Y_high*nX + iX + (1)*nX*nY*nZ) =
        *(phi_fv + iZ*nX*nY + Y_low*nX + iX + (1)*nX*nY*nZ);
}
}
/*---------------------------------------------------------------------------*/
void Periodic_XZ(PREC * phi, int nX, int nY, int nZ, int Y_to, int Y_from)
{
int iX;
int iZ;
for (iX=0; iX<nX; iX++) for (iZ=0; iZ<nZ; iZ++)
  *(phi + iZ*nX*nY + Y_to*nX + iX) =
    *(phi + iZ*nX*nY + Y_from*nX + iX);

}
/*---------------------------------------------------------------------------*/
void Periodic_XY(PREC * phi, int nX, int nY, int nZ, int Z_to, int Z_from)
{
int iX;
int iY;
for (iX=0; iX<nX; iX++) for (iY=0; iY<nY; iY++)
  *(phi + Z_to*nX*nY + iY*nX + iX) =
    *(phi + Z_from*nX*nY + iY*nX + iX);

}
/*---------------------------------------------------------------------------*/
void Periodic_YZ_flux(PREC * phi_fv, PREC * u_buff, int nX, int nY, int nZ,
  int X_low, int X_high)
{
int iY;
int iZ;
PREC u_face;

for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++) {
    u_face = LT(0.5)*( *(u_buff + X_low + iY*nX + iZ*nX*nY + (0)*nX*nY*nZ) 
      + *(u_buff + (X_low+1) + iY*nX + iZ*nX*nY + (0)*nX*nY*nZ) );
    if (u_face > 0)
      *(phi_fv + iZ*nX*nY + iY*nX + X_low + (0)*nX*nY*nZ) =
        *(phi_fv + iZ*nX*nY + iY*nX + X_high + (0)*nX*nY*nZ);
    else
      *(phi_fv + iZ*nX*nY + iY*nX + X_high + (0)*nX*nY*nZ) =
        *(phi_fv + iZ*nX*nY + iY*nX + X_low + (0)*nX*nY*nZ);
}
}
/*---------------------------------------------------------------------------*/
void Periodic_YZ(PREC * phi, int nX, int nY, int nZ, int X_to, int X_from)
{
int iY;
int iZ;
for (iY=0; iY<nY; iY++) for (iZ=0; iZ<nZ; iZ++)
  *(phi + iZ*nX*nY + iY*nX + X_to) =
    *(phi + iZ*nX*nY + iY*nX + X_from);

}
/*---------------------------------------------------------------------------*/
void Phi_3D_BC_2 (PREC * phi, int nX, int nY, int nZ, PREC delta_x)
{
Neumann_XZ(phi, nX, nY, nZ, 1,       1, LT(0.), delta_x);
Neumann_XZ(phi, nX, nY, nZ, (nY-2), -1, LT(0.), delta_x);
Neumann_XY(phi, nX, nY, nZ, 1,       1, LT(0.), delta_x);
Neumann_XY(phi, nX, nY, nZ, (nZ-2), -1, LT(0.), delta_x);
Neumann_YZ(phi, nX, nY, nZ, 1,       1, LT(0.), delta_x);
Neumann_YZ(phi, nX, nY, nZ, (nX-2), -1, LT(0.), delta_x);

Neumann_XZ(phi, nX, nY, nZ, 0,       1, LT(0.), delta_x);
Neumann_XZ(phi, nX, nY, nZ, (nY-1), -1, LT(0.), delta_x);
Neumann_XY(phi, nX, nY, nZ, 0,       1, LT(0.), delta_x);
Neumann_XY(phi, nX, nY, nZ, (nZ-1), -1, LT(0.), delta_x);
Neumann_YZ(phi, nX, nY, nZ, 0,       1, LT(0.), delta_x);
Neumann_YZ(phi, nX, nY, nZ, (nX-1), -1, LT(0.), delta_x);
}
/*---------------------------------------------------------------------------*/
void Phi_2D_BC_1 (PREC * phi, int nX, int nY, PREC delta_x)
{
Periodic_XZ(phi, nX, nY, 1, 1,      (nY-3));
Periodic_XZ(phi, nX, nY, 1, (nY-2), 2);
Periodic_YZ(phi, nX, nY, 1, 1,      (nX-3));
Periodic_YZ(phi, nX, nY, 1, (nX-2), 2);

Periodic_XZ(phi, nX, nY, 1, 0,      (nY-4));
Periodic_XZ(phi, nX, nY, 1, (nY-1), 3);
Periodic_YZ(phi, nX, nY, 1, 0,      (nX-4));
Periodic_YZ(phi, nX, nY, 1, (nX-1), 3);

}
/*---------------------------------------------------------------------------*/
void Phi_2D_BC_2 (PREC * phi, int nX, int nY, PREC delta_x)
{
Neumann_XZ(phi, nX, nY, 1, 1,       1, LT(0.), delta_x);
Neumann_XZ(phi, nX, nY, 1, (nY-2), -1, LT(0.), delta_x);
Neumann_YZ(phi, nX, nY, 1, 1,       1, LT(0.), delta_x);
Neumann_YZ(phi, nX, nY, 1, (nX-2), -1, LT(0.), delta_x);

Neumann_XZ(phi, nX, nY, 1, 0,       1, LT(0.), delta_x);
Neumann_XZ(phi, nX, nY, 1, (nY-1), -1, LT(0.), delta_x);
Neumann_YZ(phi, nX, nY, 1, 0,       1, LT(0.), delta_x);
Neumann_YZ(phi, nX, nY, 1, (nX-1), -1, LT(0.), delta_x);
}
/*---------------------------------------------------------------------------*/
void flux_2D_BC_1 (PREC * phi_fv, PREC * u_buff, int nX, int nY,
  PREC delta_x)
{
Periodic_XZ_flux(phi_fv, u_buff, nX, nY, 1, 1, (nY-3));
Periodic_XZ_flux(phi_fv, u_buff, nX, nY, 1, (nY-3), 1);
Periodic_YZ_flux(phi_fv, u_buff, nX, nY, 1, 1, (nX-3));
Periodic_YZ_flux(phi_fv, u_buff, nX, nY, 1, (nX-3), 1);
}
/*---------------------------------------------------------------------------*/
void facevals_2D_BC_1 (PREC * phi_fv, int nX, int nY,
  PREC delta_x)
{
Periodic_XZ(phi_fv+nX*nY, nX, nY, 1, 1, (nY-3));
Periodic_YZ(phi_fv      , nX, nY, 1, 1, (nX-3));
}
/*---------------------------------------------------------------------------*/
void facevals_2D_BC_2 (PREC * phi, PREC * phi_fv, int nX, int nY,
  PREC delta_x)
{
Zero_XZ(phi_fv+nX*nY, nX, nY, 1, 1);
Zero_XZ(phi_fv+nX*nY, nX, nY, 1, (nY-3));
Zero_YZ(phi_fv      , nX, nY, 1, 1);
Zero_YZ(phi_fv      , nX, nY, 1, (nX-3));
}
/*---------------------------------------------------------------------------*/

