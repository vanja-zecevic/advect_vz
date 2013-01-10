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

#ifndef FLUX_H
    #define FLUX_H

    #pragma vzg split SCHEME [up laxwend quickest_a quickest_b upmulti utopia \
      utopia_simp]

    void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha);

    #pragma vzg splitend

#endif

