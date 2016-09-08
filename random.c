/* pseudo-random number generation routines */

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
#include "random.h"

dubfloat_t rejection_sampling(){
  dubfloat_t t;
  dubfloat_t u;

  t=50.0*(dubfloat_t)rand()/RAND_MAX;
  u=(dubfloat_t)rand()/RAND_MAX;

  while(u>(dubfloat_t)exp(-((dubfloat_t)pow(t-18.0,2))/2.0)){
      t=50.0*(dubfloat_t)rand()/RAND_MAX;
      u=(dubfloat_t)rand()/RAND_MAX;
  }

  if((dubfloat_t)rand()/RAND_MAX>0.5)
    t*=-1.0;

  return t;
}
