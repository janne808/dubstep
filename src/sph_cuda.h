/* smoothed particle hydrodynamics cuda header */

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

#ifndef __SPHCUDAH__
#define __SPHCUDAH__

#ifdef __cplusplus
extern "C" void compute_smoothing_length_neighbours_cuda(struct universe *world, int iterations, int N_target);
#else
void compute_smoothing_length_neighbours_cuda(struct universe *world, int iterations, int N_target);
#endif

void smoothing_length_iterator(dubfloat_t *r, dubfloat_t *origin, int *buffer, int *buffer_index, dubfloat_t *h, int n);

#endif
