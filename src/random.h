/* pseudo-random number generation headers */

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

#ifndef __RANDOMH__
#define __RANDOMH__

static inline dubfloat_t boxmuller(){
  return sqrt(-2*log((dubfloat_t)rand()/RAND_MAX))*cos(2.0*PI*(dubfloat_t)rand()/RAND_MAX);
}

static inline dubfloat_t uniform_rand(){
  return 2.0*(dubfloat_t)rand()/RAND_MAX-1.0;
}

dubfloat_t rejection_sampling();

#endif
