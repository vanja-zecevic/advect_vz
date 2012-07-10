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
#include "src/adv/tools.h"

/*----------------------------------------------------------------------------*/
void Save_SC_XY_Slice_3D(PREC * phi, int nX, int nY, int nZ, char * out_fname)
{
/* counters */
int iX;
int iY;

FILE * out_f;

out_f = fopen(out_fname, "w");

/* Write header information */
fprintf(out_f,
  "# vtk DataFile Version 2.0\n"
  "xxx\n"
  "ASCII\n"
  "DATASET STRUCTURED_POINTS\n"
  "DIMENSIONS %i %i 1\n"
  "ORIGIN 0 0 0\n"
  "SPACING 1 1 1\n"
  "POINT_DATA %i\n"
  "SCALARS c double 1\n"
  "LOOKUP_TABLE default\n", 
  (nX), (nY), (nX)*(nY) );

for (iY=0; iY<(nY); iY++) for (iX=0; iX<(nX); iX++) {
    fprintf(out_f, "%12.11g\n", *(phi + (nZ/2)*nX*nY + iY*nX + iX));
}
fclose(out_f);
}
/*----------------------------------------------------------------------------*/
void Save_U_XY_Slice_3D(PREC * u_buff, int nX, int nY, int nZ, char * out_fname)
{
/* counters */
int iX;
int iY;
int dir = 0;

FILE * out_f;

out_f = fopen(out_fname, "w");

/* Write header information */
fprintf(out_f,
  "# vtk DataFile Version 2.0\n"
  "xxx\n"
  "ASCII\n"
  "DATASET STRUCTURED_POINTS\n"
  "DIMENSIONS %i %i 1\n"
  "ORIGIN 0 0 0\n"
  "SPACING 1 1 1\n"
  "POINT_DATA %i\n"
  "SCALARS c double 1\n"
  "LOOKUP_TABLE default\n", 
  (nX), (nY), (nX)*(nY) );

for (iY=0; iY<(nY); iY++) for (iX=0; iX<(nX); iX++) {
    fprintf(out_f, "%12.11g\n",
      *(u_buff + dir*nX*nY*nZ + (nZ/2)*nX*nY + iY*nX + iX));
}
fclose(out_f);
}
/*----------------------------------------------------------------------------*/

