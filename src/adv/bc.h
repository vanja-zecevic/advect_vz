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

#ifndef BC_H
    #define BC_H

    void Neumann_XZ(PREC * phi, int nX, int nY, int nZ, int Y_loc, int normal,
      PREC grad, PREC delta_x);

    void Neumann_XY(PREC * phi, int nX, int nY, int nZ, int Z_loc, int normal,
      PREC grad, PREC delta_x);

    void Neumann_YZ(PREC * phi, int nX, int nY, int nZ, int X_loc, int normal,
      PREC grad, PREC delta_x);

    void Zero_XZ(PREC * phi, int nX, int nY, int nZ, int Y_loc);
    void Zero_XY(PREC * phi, int nX, int nY, int nZ, int Z_loc);
    void Zero_YZ(PREC * phi, int nX, int nY, int nZ, int X_loc);

    void Periodic_XZ_fv(PREC * phi_fv, PREC * u_buff, int nX, int nY, int nZ,
      int Y_low, int Y_high);
    void Periodic_XZ(PREC * phi, int nX, int nY, int nZ, int Y_to, int Y_from);
    void Periodic_XY(PREC * phi, int nX, int nY, int nZ, int Z_to, int Z_from);
    void Periodic_YZ_fv(PREC * phi_fv, PREC * u_buff, int nX, int nY, int nZ,
      int X_low, int X_high);
    void Periodic_YZ(PREC * phi, int nX, int nY, int nZ, int X_to, int X_from);

    void Phi_3D_BC_2 (PREC * phi, int nX, int nY, int nZ, PREC delta_x);
    void Phi_2D_BC_1 (PREC * phi, int nX, int nY, PREC delta_x);
    void Phi_2D_BC_2 (PREC * phi, int nX, int nY, PREC delta_x);
    void flux_2D_BC_1 (PREC * phi_fv, PREC * u_buff, int nX, int nY,
      PREC delta_x);
    void facevals_2D_BC_1 (PREC * phi_fv, int nX, int nY,
      PREC delta_x);
    void facevals_2D_BC_2 (PREC * phi, PREC * phi_fv, int nX, int nY,
      PREC delta_x);

#endif
