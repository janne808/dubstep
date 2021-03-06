
/* barnes-hut tree headers */

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

#ifndef __TREEH__
#define __TREEH__

/* tree cell */
struct cell{ 
  dubfloat_t *r; /* cell particle displacement vector */
  dubfloat_t *m; /* cell particle mass vector */

  int *particle_index_list; /* cell particle index relative to the world vector */

  int particle_index; /* particle index for tree branch end cell */

  dubfloat_t space[3*8]; /* cell spatial limits */
    
  dubfloat_t l;   /* cell size */
    
  dubfloat_t mass; /* total mass of particles */
  int num; /* number of particles in cell */
  dubfloat_t center[3]; /* cell center of mass */
  dubfloat_t distr_len; /* mass distribution length */
    
  int parent; /* parent cell index */
  int children[8]; /* children cell indices */
  int numcells; /* amount of cells in tree */
  int index;
  int numchild;
};

/* function prototypes */
int init_treeroot(struct cell *tree, struct universe *world, dubfloat_t *r);

void direct_summation(struct universe *world, dubfloat_t *r, dubfloat_t *a, dubfloat_t G);

void force_walk(struct universe *world, struct cell *tree, struct cell *root, dubfloat_t *r,
		dubfloat_t *a, dubfloat_t G, dubfloat_t theta, dubfloat_t h);

void neighbour_walk(struct cell *tree, struct cell *root, dubfloat_t *r, dubfloat_t h, dubfloat_t max_h,
		    dubfloat_t *h_in, int *neighbour_num, int *neighbour_list);

void compute_total_energy(struct universe *world, dubfloat_t theta, int lo, int hi);

void tree_branch(struct universe *world, struct cell *tree, struct cell *root, int *cellindex);

void tree_recurse(struct cell *tree, struct cell *root);

void neighbour_recurse(struct cell *tree, struct cell *root, dubfloat_t *r, dubfloat_t h, dubfloat_t max_h,
		       dubfloat_t *h_in, int *neighbour_num, int *neighbour_list);

void potential_recurse(struct cell *tree, struct cell *root, 
		       dubfloat_t *r, dubfloat_t m, dubfloat_t *U, dubfloat_t G, dubfloat_t theta,
		       dubfloat_t epsilon);

void force_recurse(struct cell *tree, struct cell *root,
		   dubfloat_t *r, dubfloat_t *f, dubfloat_t G, dubfloat_t theta,
		   dubfloat_t epsilon);

void branch_recurse(struct universe *world, struct cell *tree, struct cell *root, int *cellindex);

#endif
