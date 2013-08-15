/* dubstep headers */

/*
 *  (C) 2013 Janne Heikkarainen <janne.heikkarainen@tut.fi>
 *
 *  All rights reserved.
 *
 *  This file is part of Dubstep POSIX/CUDA Self-gravitating Smoothed Particle Hydrodynamics Simulator.
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
#define ENABLE_GUI 0
#endif

#ifndef ADAPTIVE_SMOOTHING
#define ADAPTIVE_SMOOTHING 0
#endif

#ifndef MIN_SMOOTH_LEN
#define MIN_SMOOTH_LEN 0.7
#endif

#ifndef MAX_SMOOTH_LEN
#define MAX_SMOOTH_LEN 0.7
#endif

#ifndef ADAPTIVE_SOFTENING
#define ADAPTIVE_SOFTENING 0
#endif

#ifndef SOFTENING_FACTOR
#define SOFTENING_FACTOR 0.7
#endif

#ifndef GEOMETRIC_MEAN_SYMMETRIZATION
#define GEOMETRIC_MEAN_SYMMETRIZATION 0
#endif

#ifndef ARITHMETIC_MEAN_SYMMETRIZATION
#define ARITHMETIC_MEAN_SYMMETRIZATION 1
#endif

#ifndef NUM
#define NUM 10000
#endif

#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif

#ifndef PI
#define PI 3.14159265358979323846264
#endif

#if CUDA
typedef float dubfloat_t;
#else
typedef double dubfloat_t;
#endif

/* structure for RGB color */
struct color{
  dubfloat_t r;
  dubfloat_t g;
  dubfloat_t b;
};

/* structure for list of interacting particles */
struct particlelist{
  int num;
  int max_size;
  int *list;
};

/* universe state structure */
struct universe{
  int dim; /* number of dimensions */
  int num; /* number of particles */

  int *num_neighbours; /* number of neighbouring particles */
  struct particlelist *neighbour_list; /* pointer list of list of interacting particles */

  struct cell *tree; /* tree cell array */

  dubfloat_t G; /* gravitational constant */

  dubfloat_t *r; /* particle displacement vectors */
  dubfloat_t *r2; /* particle displacement halfstep vectors */
  dubfloat_t *v; /* particle velocity vectors */ 
  dubfloat_t *v2; /* particle velocity halfstep vectors */ 
  dubfloat_t *a; /* particle acceleration vectors*/
  dubfloat_t *a2; /* particle acceleration halfstep vectors*/

  dubfloat_t *m; /* particle mass array */

  dubfloat_t *rho; /* particle density array */
  dubfloat_t *p; /* particle pressure array */
  dubfloat_t *c; /* particle speed of sound array */

  dubfloat_t *u; /* particle internal energy array */
  dubfloat_t *u2; /* particle internal energy halfstep array */
  dubfloat_t *du; /* particle internal energy array */
  dubfloat_t *du2; /* particle internal energy halfstep array */

  dubfloat_t *del_rho; /* particle density partial derivative array */

  dubfloat_t *h; /* particle smoothing length parameter array */

  dubfloat_t *dt_CFL; /* particle time step array */

  int *cellindex; /* particle's residing tree cell index */

  dubfloat_t *a_sph; /* particle hydrodynamic acceleration vector array */
  dubfloat_t *a_tree; /* particle gravitational acceleration vector array */

  /* model parameters */
  dubfloat_t dt;
  dubfloat_t epsilon;
  dubfloat_t gamma;
  dubfloat_t theta;
  dubfloat_t alpha;
  dubfloat_t beta;

  dubfloat_t time; /* current universal time */

  int timediv;
  dubfloat_t sub_dt;

  int *kick;
  int *kick_list;
  int kick_num;

  int *time_bin;

  /* total system energy */
  dubfloat_t u_int;
  dubfloat_t u_grav;
  dubfloat_t u_kin;

  dubfloat_t avg_u;
  dubfloat_t min_u;
};

#endif
