/* smoothed particle hydrodynamics header */

/*
 *  (C) 2012 Janne Heikkarainen <janne.heikkarainen@tut.fi>
 *
 *  All rights reserved.
 *
 *  This file is part of Dubstep ANSI C Self-gravitating Smoothed Particle Hydrodynamics Simulator.
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

double kernel(double q, double h);
double kernel_d(double q, double h);

double artificial_viscosity(double *dv_ij, double h_ij, double rho_ij,
			    double c_ij, double *dr, double rr, double alpha, double beta, double neta);

void compute_smoothing_length_neighbours(struct universe *world, int iterations, int N_target, int lo, int hi);
void compute_smoothing_length_mass(struct universe *world, int iterations, double m_target, int lo, int hi);
void compute_density(struct universe *world, double h, int lo, int hi);
void compute_pressure(struct universe *world, double gamma, int lo, int hi);
void compute_soundspeed(struct universe *world, double gamma, int lo, int hi);
void compute_cfl(struct universe *world, double C_0, int lo, int hi);
void compute_internal_energy(struct universe *world, double *a_sph, int lo, int hi);
void compute_sph_acceleration(struct universe *world, double *r, double *v, double *a, int lo, int hi);

#endif
