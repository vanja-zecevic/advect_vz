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

#ifndef INIT_H
    #define INIT_H

    void init_u_rigid(PREC * u, int nX, int nY, int nZ, PREC u_init_max);

    void init_u_const(PREC * u, int nX, int nY, int nZ, PREC u_init_max,
      PREC u_init_angle);

    void init_u_tvort(PREC * u_buff, int nX, int nY, int nZ, PREC u_init_max);

    void init_phi_cyl(PREC * phi, int nX, int nY, int nZ, PREC x_zero,
      PREC y_zero, PREC z_zero, PREC phi_init_max, PREC radius);

    void init_phi_gauss(PREC * phi, int nX, int nY, int nZ, PREC x_zero,
      PREC y_zero, PREC z_zero, PREC phi_init_max, PREC B);

    void init_phi_step(PREC * phi, int nX, int nY, int nZ, PREC phi_init_max,
      PREC h_step);

    void init_phi_const(PREC * phi, int nX, int nY,
      int nZ, PREC phi_init_max);

#endif
