/* linear algebra headers */

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

#ifndef __LINEARH__
#define __LINEARH__

static inline dubfloat_t euclidean_norm(dubfloat_t *r, int dim){
  /* loop variable */
  int ii;

  /* inner product sum variable */
  dubfloat_t sum;

  sum=0;
  for(ii=0;ii<dim;ii++){
    sum+=r[ii]*r[ii];
  }

  return sqrt(sum);
}

static inline dubfloat_t euclidean_distance(dubfloat_t *r, dubfloat_t *v, int dim){
  /* loop variable */
  int ii;

  /* inner product sum variable */
  dubfloat_t sum;

  sum=0;
  for(ii=0;ii<dim;ii++){
    sum+=(r[ii]-v[ii])*(r[ii]-v[ii]);
  }

  return sqrt(sum);
}

#endif
