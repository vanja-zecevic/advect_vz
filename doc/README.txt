Copyright (C) 2012 Vanja Zecevic
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.

================================================================================
Introduction:
                              --------------------------------------------------
This program was designed in order to test various implementations of finite
volume difference schemes and flux limiters. The code simulates the advection
and diffusion of a scalar with a specified, velocity field. Several test cases
including the advection of a Gaussian hill and a deformational flow are
provided.

================================================================================
Author:
                              --------------------------------------------------
This software was written by, 

Vanja Zecevic
School of Aerospace, Mechanical and Mechatronic Engineering
The University of Sydney
Sydney, NSW, 2006
Australia

Please send any correspondence to,

vanja.zecevic@aeromech.usyd.edu.au

================================================================================
Use:
                              --------------------------------------------------
First compile the program

cd advect_vz
make

then create a run directory with a config file and the binary,

mkdir run
cp doc/project.conf.hill ./run/project.conf
cp bin/run2d ./run
cd run
./run2d

You can change various parameters including the type of simulation to run in the
configuration file. Users are encouraged to edit the source files for further
experimentation.

The program outputs .vtk files of the scalar field at a specified interval. You
will need a scientific visualization program such as visit
  https://wci.llnl.gov/codes/visit/
or paraview
  http://www.paraview.org/
in order to view the data.

================================================================================
Difference schemes and algorithms:
                              --------------------------------------------------

This program uses a code generator (gen_vz) in order to avoid branching code in 
the lowest level functions. Several alternative algorithms are also included to
illustrate the advantage of using the code generator. From the config file,

# algo
# 0 = generated files.                 1 = branching, lowest level.
# 2 = branching, higher level.         3 = function pointers.
# 4 = lliffe vectors, only quickest.   5 = face values
# 6 = multi-dimensional flux limiter

Algorithms 0, 1, 2 and 3 are similar apart from using different methods for
branching at point where face values are calculated. Algorithm 4 uses lliffe
vectors instead of strided access. Algorithm 5 updates the cell in two passes,
the first pass calculates and writes each face value, the second pass updates
the cell values from stored face values. This prevents calculating face values
twice. Algorithm 6 uses the multi-dimensional flux limiter [4].

The following difference schemes are included.

# scheme
# 0 = upwind         1 = laxwend
# 2 = quickest_a     3 = quickest_b
# 4 = upmulti        5 = utopia_simp
# 6 = utopia

There are two implementations of the quickest scheme, quickest_a knows all 4
memory access locations at compile time while quickest_b dynamically calculates
the offsets required to access only the 3 required cells. The utopia_simp
algorithm is a modification of utopia without 'twist' terms, it represents a
flux integral using second order interpolation in the upstream region of the
face. The same interpolation formula is used for the entire upwind parallelogram
regardless of which cells it intersects, in contrast, utopia performs a separate
interpolation for each upwind cell. 

================================================================================
Reference material:
                              --------------------------------------------------
Papers describing the difference schemes and flux limiter used in this program.

[1] Leonard, B. P.
Stable and accurate convective modeling procedure based on quadratic upstream
  interpolation
Comput. Meth. Appl. Mech. Eng. 19 (1979)

[2]  Leonard, B. P.
The ultimate conservative difference scheme applied to unsteady one-dimensional
  advection
Comput. Meth. Appl. Mech. Eng. 88 (1991)

[3] Leonard, B. P. and MacVean, M. K. and Lock, A. P.
The flux integral method for multidimensional convection and diffusion
Appl. Math. Model. 19 (1995)

[4] Thuburn, J
Multidimensional flux-limited advection schemes
J. Comput. Phys. 123 (1996)

