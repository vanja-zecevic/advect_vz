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

#ifndef CORE_TEST_H
    #define CORE_TEST_H

    void update_2D_quickest_a_lliffe (PREC ** phi_a, PREC ** phi_b, 
      PREC * u_buff, int iT, int nX,
      int nY, int bc_phi, PREC delta_t, PREC alpha);

    inline void update_cell_2D_quickest_a_lliffe (PREC ** sc_old, PREC ** sc_new,
      PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha);

    void update_2D_BRANCH(PREC * phi_a, PREC * phi_b, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha, int scheme);

    void update_cell_2D_BRANCH(PREC * sc_old, PREC * sc_new, PREC * u_buff,
      int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha, int scheme);

    void update_2D_BRANCH2(PREC * phi_a, PREC * phi_b, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha, int scheme);

    void update_cell_2D_BRANCH2(PREC * sc_old, PREC * sc_new, PREC * u_buff,
      int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha, int scheme);

    void update_2D_PTRS(PREC * phi_a, PREC * phi_b, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha, int scheme);

    void update_cell_2D_PTRS(PREC * sc_old, PREC * sc_new, PREC * u_buff,
      int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha, int scheme,
      void (*get_face_2d_ptr) (int iX, int iY, int nX, int nY, int X_nbr,
      int Y_nbr, int nrm, int tid, PREC delta_t, PREC u_face, PREC dummy,
      PREC * sc_face, PREC * sc_old));

#endif

