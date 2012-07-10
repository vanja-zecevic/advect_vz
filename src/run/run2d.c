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
#include "src/adv/bc.h"
#include "src/adv/flux.h"
#include "src/adv/core.h"
#include "src/adv/core_test.h"
#include "src/adv/tools.h"
#include "libvz/exit_vz.h"
#include "libvz/cfg_vz.h"
#include "libvz/lliffe_vz.h"

int main(int argc, char *argv[])
{
int iT;
char out_fname[256];
char * projfilename = "project.conf";
struct timespec t_start;
struct timespec t_end;
struct timespec t_d1;
struct timespec t_d2;
PREC MUPS;
PREC * u_buff = NULL;

lliffe_2d_(PREC) phi_lliffe_a_struct;
PREC ** phi_lliffe_a = NULL;
PREC * phi_a = NULL;

lliffe_2d_(PREC) phi_lliffe_b_struct;
PREC ** phi_lliffe_b = NULL;

PREC * phi_b = NULL;
PREC * phi_fv = NULL;

void (*update_2D_nrm)  (PREC * phi_a, PREC * phi_b, PREC * u_buff, int iT,
  int nX, int nY, int bc_phi, PREC delta_t, PREC alpha) = NULL;
void (*update_2D_face) (PREC * phi_a, PREC * phi_fv, PREC * u_buff, int iT,
  int nX, int nY, int bc_phi, PREC delta_t, PREC alpha) = NULL;
void (*update_2D_flux) (PREC * phi_a, PREC * phi_fv, PREC * u_buff, int iT,
  int nX, int nY, int bc_phi, PREC delta_t, PREC alpha) = NULL;

#include "libvz/cfg_def.inc" /* Config file variables.  */

#include "libvz/cfg_parse.inc" /* Parse the config file.  */

if      (scheme==0) {
    update_2D_nrm  = &update_2D_nrm_up;
    update_2D_face = &update_2D_face_up;
    update_2D_flux = &update_2D_flux_up;
}
else if (scheme==1) {
    update_2D_nrm  = &update_2D_nrm_laxwend;
    update_2D_face = &update_2D_face_laxwend;
    update_2D_flux = &update_2D_flux_laxwend;
}
else if (scheme==2) {
    update_2D_nrm  = &update_2D_nrm_quickest_a;
    update_2D_face = &update_2D_face_quickest_a;
    update_2D_flux = &update_2D_flux_quickest_a;
}
else if (scheme==3) {
    update_2D_nrm  = &update_2D_nrm_quickest_b;
    update_2D_face = &update_2D_face_quickest_b;
    update_2D_flux = &update_2D_flux_quickest_b;
}
else if (scheme==4) {
    update_2D_nrm  = &update_2D_nrm_upmulti;
    update_2D_face = &update_2D_face_upmulti;
    update_2D_flux = &update_2D_flux_upmulti;
}
else if (scheme==5) {
    update_2D_nrm  = &update_2D_nrm_utopia_simp;
    update_2D_face = &update_2D_face_utopia_simp;
    update_2D_flux = &update_2D_flux_utopia_simp;
}
else if (scheme==6) {
    update_2D_nrm  = &update_2D_nrm_utopia;
    update_2D_face = &update_2D_face_utopia;
    update_2D_flux = &update_2D_flux_utopia;
}
else
    exit_msg("Invalid scheme.", 1);

/* Allocate u and phi.  */
u_buff = (PREC*)malloc(3*nX*nY*sizeof(PREC));
if (u_buff==NULL) return 1;

if (algo==0 || algo==1 || algo==2 || algo==3) {
    phi_a = (PREC*)malloc(nY*nX*sizeof(PREC));
    phi_b = (PREC*)malloc(nY*nX*sizeof(PREC));
}
if (algo==4) {
    malloc_2d_(PREC) (&phi_lliffe_a_struct, nY, nX);
    phi_lliffe_a = phi_lliffe_a_struct.sec;
    phi_a = phi_lliffe_a_struct.pri;

    malloc_2d_(PREC) (&phi_lliffe_b_struct, nY, nX);
    phi_lliffe_b = phi_lliffe_b_struct.sec;
}
if (algo==5 || algo==6) {
    phi_a = (PREC*)malloc(nY*nX*sizeof(PREC));
    phi_fv = (PREC*)malloc(2*nY*nX*sizeof(PREC));
}
/* Initialise u and phi.  */
if      (u_init_type == 0)
  init_u_rigid(u_buff, nX, nY, 1, u_init_max);
else if (u_init_type == 1)
  init_u_const(u_buff, nX, nY, 1, u_init_max, u_init_const_angle);
else if (u_init_type == 2)
  init_u_tvort(u_buff, nX, nY, 1, u_init_max);

if      (phi_init_type == 0)
  init_phi_const(phi_a, nX, nY, 1, phi_init_max);
if      (phi_init_type == 1)
  init_phi_gauss(phi_a, nX, nY, 1, ((PREC)(nX-1)/LT(2.)), gauss_offset,
    ((PREC)(nZ-1)/LT(2.)), phi_init_max, gauss_B);
else if (phi_init_type == 2)
  init_phi_step (phi_a, nX, nY, 1, phi_init_max, phi_init_h_step);
else if (phi_init_type == 3)
  init_phi_cyl  (phi_a, nX, nY, 1, ((PREC)(nX-1)/LT(2.)), gauss_offset,
    ((PREC)(nZ-1)/LT(2.)), phi_init_max, phi_init_cyl_rad);

if      (bc_phi == 1)
    Phi_2D_BC_1(phi_a, nX, nY, DELTA_X);
else if (bc_phi == 2)
    Phi_2D_BC_2(phi_a, nX, nY, DELTA_X);

/* Sets the number of parallel openMP threads to use */
omp_set_num_threads(num_threads);

printf("Solver starting...\n");
clock_gettime(CLOCK_REALTIME, &t_start);
t_d1.tv_sec = t_start.tv_sec;
t_d1.tv_nsec = t_start.tv_nsec;

for (iT=0; iT<nT; iT++) {

    if      (algo==0)
        update_2D_nrm(phi_a, phi_b, u_buff,
          iT, nX, nY, bc_phi, delta_t, alpha);
    else if (algo==1)
        update_2D_BRANCH(phi_a, phi_b, u_buff,
          iT, nX, nY, bc_phi, delta_t, alpha, scheme);
    else if (algo==2)
        update_2D_BRANCH2(phi_a, phi_b, u_buff,
          iT, nX, nY, bc_phi, delta_t, alpha, scheme);
    else if (algo==3)
        update_2D_PTRS(phi_a, phi_b, u_buff,
          iT, nX, nY, bc_phi, delta_t, alpha, scheme);
    else if (algo==4)
        update_2D_quickest_a_lliffe(phi_lliffe_a, phi_lliffe_b, u_buff,
          iT, nX, nY, bc_phi, delta_t, alpha);
    else if (algo==5)
        update_2D_face(phi_a, phi_fv, u_buff,
          iT, nX, nY, bc_phi, delta_t, alpha);
    else if (algo==6)
        update_2D_flux(phi_a, phi_fv, u_buff,
          iT, nX, nY, bc_phi, delta_t, alpha);

    if (iT%t_report==0) {
        clock_gettime(CLOCK_REALTIME, &t_d2);
        MUPS = LT(1e-6)*(float)(nX*nY*1*t_report)
         /( (long)t_d2.tv_sec - (long)t_d1.tv_sec
         +(t_d2.tv_nsec-t_d1.tv_nsec)/LT(1000000000.) );
        printf("i = %5.5i, %g MUPS\n", iT, MUPS);
        t_d1.tv_sec = t_d2.tv_sec;
        t_d1.tv_nsec = t_d2.tv_nsec;
    }
    if (iT%t_save==0) {
        snprintf(out_fname,256,"cube%i.vtk", iT);
        *(out_fname + 255) = '\0';
        Save_SC_XY_Slice_3D(phi_a, nX, nY, nZ, out_fname);
    }
}

clock_gettime(CLOCK_REALTIME, &t_end);
printf("Solver finished, time elapsed\n%g s\n",
  ((double)t_end.tv_sec-(double)t_start.tv_sec) +
  ((double)t_end.tv_nsec-(double)t_start.tv_nsec)/LT(1000000000.));

return 0;
}
/*----------------------------------------------------------------------------*/

