/* barnes-hut tree headers */

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

#ifndef __TREEH__
#define __TREEH__

/* tree cell */
struct cell{ 
  double *r; /* cell particle displacement vector */
  double *m; /* cell particle mass vector */

  int *particle_index_list; /* cell particle index relative to the world vector */

  int particle_index; /* particle index for tree branch end cell */

  double space[3*8]; /* cell spatial limits */
    
  double l;   /* cell size */
    
  double mass; /* total mass of particles */
  int num; /* number of particles in cell */
  double center[3]; /* cell center of mass */
    
  int parent; /* parent cell index */
  int children[8]; /* children cell indices */
  int numcells; /* amount of cells in tree */
  int index;
  int numchild;
};

/* function prototypes */
int init_treeroot(struct cell *tree, struct universe *world, double *r);
double compute_total_potential_energy(struct universe *world);
void treebranch(struct cell *tree, struct cell *root, int *cellindex);
void treerecurse(struct cell *tree, struct cell *root);
void neighbourrecurse(struct cell *tree, struct cell *root, double *r, double h, double max_h,
		      double *h_in, int *neighbour_num, int *neighbour_list);
void potentialrecurse(struct cell *tree, struct cell *root, 
		      double *r, double m, double *U, double G, double theta,
		      double epsilon);
void forcerecurse(struct cell *tree, struct cell *root,
		  double *r, double *f, double G, double theta,
		  double epsilon);
void branchrecurse(struct cell *tree, struct cell *root, int *cellindex);

#endif
