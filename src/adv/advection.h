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

#ifndef ADVECTION_H
    #define ADVECTION_H

    #include <math.h>
    #include <time.h>
    #include <omp.h>
    #include <stdio.h>
    #include <stdlib.h>

    #define PREC double
    #define LT(literal) literal

    #define DELTA_X LT(1.)
    #define DELTA_X_INV LT(1.)
    #define VZ_PI 3.14159265358979323846
    #define THIRD 0.33333333333333333333
    #define CVARS_FILE "src/adv/cvars.def"

#endif
