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

    #define SCHEME up
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_up
#define SCHEME_NUM 0

    void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha);

    inline void APND_SCHEME(update_fvflux_2D) (PREC * sc_old, PREC * sc_fv,
      PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha);

    
#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME laxwend
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_laxwend
#define SCHEME_NUM 1

    void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha);

    inline void APND_SCHEME(update_fvflux_2D) (PREC * sc_old, PREC * sc_fv,
      PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha);

    
#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME quickest_a
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_quickest_a
#define SCHEME_NUM 2

    void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha);

    inline void APND_SCHEME(update_fvflux_2D) (PREC * sc_old, PREC * sc_fv,
      PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha);

    
#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME quickest_b
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_quickest_b
#define SCHEME_NUM 3

    void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha);

    inline void APND_SCHEME(update_fvflux_2D) (PREC * sc_old, PREC * sc_fv,
      PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha);

    
#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME upmulti
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_upmulti
#define SCHEME_NUM 4

    void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha);

    inline void APND_SCHEME(update_fvflux_2D) (PREC * sc_old, PREC * sc_fv,
      PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha);

    
#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME utopia
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_utopia
#define SCHEME_NUM 5

    void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha);

    inline void APND_SCHEME(update_fvflux_2D) (PREC * sc_old, PREC * sc_fv,
      PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha);

    
#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM

#define SCHEME utopia_simp
#define APND_SCHEME(IN_TOKEN) IN_TOKEN##_utopia_simp
#define SCHEME_NUM 6

    void APND_SCHEME(update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff,
      int iT, int nX, int nY, int bc_phi, PREC delta_t, PREC alpha);

    inline void APND_SCHEME(update_fvflux_2D) (PREC * sc_old, PREC * sc_fv,
      PREC * u_buff, int iX, int iY, int nX, int nY, PREC delta_t, PREC alpha);

    
#undef SCHEME
#undef APND_SCHEME
#undef SCHEME_NUM



#endif

