/* dubstep headers */

/*
 *  (C) 2012 Janne Heikkarainen <janne.heikkarainen@tut.fi>
 *
 *  All rights reserved.
 *
 *  This file is part of Dubstep ANSI C/CUDA Self-gravitating Smoothed Particle Hydrodynamics Simulator.
 *
 *  Dubstep is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Dubstep is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Dubstep.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __DUBSTEPH__
#define __DUBSTEPH__

#ifndef ENABLE_GUI
#define ENABLE_GUI 1
#endif

#define NUM 800
#define NUM_THREADS 2

#ifndef PI
#define PI 3.14159265358979323846264
#endif

/* structure for RGB color */
struct color{
  double r;
  double g;
  double b;
};

/* structure for list of interacting particles */
struct particlelist{
  int num;
  int *list;
};

/* universe state structure */
struct universe{
  int dim; /* number of dimensions */
  int num; /* number of particles */

  int *num_neighbours; /* number of neighbouring particles */
  struct particlelist *neighbour_list; /* pointer list of list of interacting particles */

  struct cell *tree; /* tree cell array */

  double G; /* gravitational constant */

  double *r; /* particle displacement vectors */
  double *r2; /* particle displacement halfstep vectors */
  double *v; /* particle velocity vectors */ 
  double *v2; /* particle velocity halfstep vectors */ 
  double *a; /* particle acceleration vectors*/
  double *a2; /* particle acceleration halfstep vectors*/

  double *m; /* particle mass array */

  double *rho; /* particle density array */
  double *p; /* particle pressure array */
  double *c; /* particle speed of sound array */

  double *u; /* particle internal energy array */
  double *u2; /* particle internal energy halfstep array */
  double *du; /* particle internal energy array */
  double *du2; /* particle internal energy halfstep array */

  double *del_rho; /* particle density partial derivative array */

  double *h; /* particle smoothing length parameter array */

  double *dt_CFL; /* particle time step array */

  int *cellindex; /* particle's residing tree cell index */

  double *a_sph; /* particle hydrodynamic acceleration vector array */
  double *a_tree; /* particle gravitational acceleration vector array */

  /* model parameters */
  double dt;
  double epsilon;
  double gamma;
  double theta;
  double alpha;
  double beta;

  double time; /* current universal time */

  int timediv;
  double sub_dt;

  int *kick;

  int *time_bin;

  /* total system energy */
  double u_int;
  double u_grav;
  double u_kin;

  double avg_u;
  double min_u;
};

#endif
