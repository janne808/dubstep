/* particle lattice routines */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dubstep.h"
#include "sph.h"
#include "tree.h"
#include "linear.h"
#include "lattice.h"

void build_particle_lattice(struct universe *world){
  /* loop variables */
  int ii;

  /* cell size */
  dubfloat_t l=1.0;

  /* maximum displacement */
  dubfloat_t max_r;

  /* number of cells */
  int cell_num;

  /* buffers */
  int *x_buffer;
  int *y_buffer;
  int *z_buffer;

  /* search for largest displacement */
  max_r=0;
  for(ii=0;ii<world->num;ii++){
    if(world->r[3*ii+0]>max_r)
      max_r=world->r[3*ii+0];
    if(world->r[3*ii+1]>max_r)
      max_r=world->r[3*ii+1];
    if(world->r[3*ii+2]>max_r)
      max_r=world->r[3*ii+2];
  }

  /* make sure lattice size is integer multiply of cell size */
  /* number of cells along one dimension */
  cell_num=1;
}

