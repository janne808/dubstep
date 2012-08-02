/* smoothed particle hydrodynamics header */

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

#ifndef __SPHH__
#define __SPHH__

void smooth_velocity_field(struct universe *world, int lo, int hi);
void smooth_energy_field(struct universe *world, int lo, int hi);

void compute_smoothing_length_tree(struct universe *world, double min_h, double max_h, int iterations,
				   int N_target, double *r, struct cell *tree, struct cell *root, int lo, int hi);
void compute_constant_smoothing_length_tree(struct universe *world, double min_h, double max_h, int iterations,
					    int N_target, double *r, struct cell *tree, struct cell *root, int lo, int hi);
void compute_smoothing_length_neighbours(struct universe *world, int iterations, int N_target, int lo, int hi);
void compute_smoothing_length_mass(struct universe *world, int iterations, double m_target, int lo, int hi);
void compute_density(struct universe *world, double h, int lo, int hi);
void compute_pressure(struct universe *world, double gamma, int lo, int hi);
void compute_soundspeed(struct universe *world, double gamma, int lo, int hi);
void compute_cfl(struct universe *world, double C_0, int lo, int hi);
void compute_internal_energy_and_acceleration(struct universe *world, double *r, double *v, double *a, int lo, int hi);
void update_time_bins(struct universe *world, int lo, int hi);

#endif
