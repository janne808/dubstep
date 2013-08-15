/* smoothed particle hydrodynamics header */

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

#ifndef __SPHH__
#define __SPHH__

static inline dubfloat_t kernel(dubfloat_t q, dubfloat_t h){
  dubfloat_t val=q;

  if(q<1.0){
    val=1.0-3.0/2.0*val*val+3.0/4.0*val*val*val;
  }
  else if(q<2.0){
    val=2.0-val;
    val=1.0/4.0*val*val*val;
  }
  else
    val=0;

  return val/(PI*h*h*h);
}

static inline dubfloat_t kernel_d(dubfloat_t q, dubfloat_t h){
  dubfloat_t val=q;

  if(q<1.0){
    val=-3.0*val+9.0/4.0*val*val;
  }
  else if(q<2.0){
    val=2.0-val;
    val=-3.0/4.0*val*val;
  }
  else
    val=0;

  return val/(PI*h*h*h);
}

static inline dubfloat_t kernel_grav_f(dubfloat_t r, dubfloat_t epsilon){
  dubfloat_t val;
  dubfloat_t u;
  dubfloat_t u2;

  u=r/epsilon;

  if(u<1.0){
    u2=u*u;
    val=-2.0/epsilon*(1.0/3.0*u2-3.0/20.0*u2*u2+1.0/20.0*u2*u2*u)+7.0/5.0*epsilon;
  }
  else if(u<2.0){
    u2=u*u;
    val=-1.0/15.0*r-1/epsilon*(4.0/3.0*u2-u2*u+3.0/10.0*u2*u2-1.0/30.0*u2*u2*u)+8.0/5.0*epsilon;
  }
  else
    val=1/r;
  
  return val;
}

static inline dubfloat_t kernel_grav_g(dubfloat_t r, dubfloat_t epsilon){
  dubfloat_t val;
  dubfloat_t u;
  dubfloat_t u2;

  u=r/epsilon;

  if(u<1.0){
    u2=u*u;
    val=1.0/(epsilon*epsilon*epsilon)*(4.0/3.0-6.0/5.0*u2+1.0/2.0*u2*u);
  }
  else if(u<2.0){
    u2=u*u;
    val=1.0/(r*r*r)*(-1.0/15.0+8.0/3.0*u2*u-3.0*u2*u2+6.0/5.0*u2*u2*u-1.0/6.0*u2*u2*u2);
  }
  else
    val=1/(r*r*r);
  
  return val;
}

static inline dubfloat_t artificial_viscosity(dubfloat_t *dv_ij, dubfloat_t h_ij, dubfloat_t rho_ij,
					  dubfloat_t c_ij, dubfloat_t *dr, dubfloat_t rr, dubfloat_t alpha, dubfloat_t beta, dubfloat_t neta){
  dubfloat_t mu_ij;

  mu_ij=0;
  if((dv_ij[0]*dr[0]+dv_ij[1]*dr[1]+dv_ij[2]*dr[2])<0.0)
    mu_ij=(dv_ij[0]*dr[0]+dv_ij[1]*dr[1]+dv_ij[2]*dr[2])/(rr*rr/(h_ij)+neta);

  return (-alpha*mu_ij*c_ij+beta*mu_ij*mu_ij)/rho_ij;
}

void build_particle_lattice(struct universe *world);

void smooth_velocity_field(struct universe *world, int lo, int hi);

void smooth_energy_field(struct universe *world, int lo, int hi);

void compute_smoothing_length_tree(struct universe *world, dubfloat_t min_h, dubfloat_t max_h, int iterations,
				   int N_target, dubfloat_t *r, struct cell *tree, struct cell *root,
				   int lo, int hi, int *buffer);

void compute_constant_smoothing_length_tree(struct universe *world, dubfloat_t min_h, dubfloat_t max_h, int iterations,
					    int N_target, dubfloat_t *r, struct cell *tree, struct cell *root,
					    int lo, int hi, int *buffer);

void compute_smoothing_length_neighbours(struct universe *world, int iterations, int N_target, int lo, int hi);

void compute_density(struct universe *world, dubfloat_t h, int lo, int hi);

void compute_pressure(struct universe *world, dubfloat_t gamma, int lo, int hi);

void compute_soundspeed(struct universe *world, dubfloat_t gamma, int lo, int hi);

void compute_cfl(struct universe *world, dubfloat_t C_0, int lo, int hi);

void compute_internal_energy_and_acceleration(struct universe *world, dubfloat_t *r, dubfloat_t *v, dubfloat_t *a, int lo, int hi);

void update_time_bins(struct universe *world, int lo, int hi);

void update_kick_list(struct universe *world);

void init_kick_list(struct universe *world);

#endif
