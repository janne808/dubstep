#ifndef __DUBSTEPH__
#define __DUBSTEPH__

#ifndef ENABLE_GUI
#define ENABLE_GUI 1
#endif

#define NUM 500
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

  double *last_kick;
};

#endif
