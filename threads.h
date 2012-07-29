/* thread headers */

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

#ifndef __THREADSH__
#define __THREADSH__

struct thread_data{
  int thread_id;
  struct universe *world;
  double var;
  int lo;
  int hi;
};

struct thread_data2{
  int thread_id;
  struct universe *world;
  double *r;
  double *v;
  double *a;
  int lo;
  int hi;
};

struct thread_data3{
  int thread_id;
  struct universe *world;
  struct cell *tree;
  double *a_sph;
  double var1;
  double var2;
  double var3;
  int lo;
  int hi;
};

struct thread_data4{
  int thread_id;
  struct universe *world;
  int var1;
  int var2;
  int lo;
  int hi;
};

struct thread_data5{
  int thread_id;
  struct universe *world;
  struct cell *tree;
  double *a_sph;
  double var1;
  double var2;
  double var3;
  int lo;
  int hi;
};

struct thread_data6{
  int thread_id;
  struct universe *world;
  struct cell *tree;
  struct cell *root;
  double *r;
  double var1;
  double var2;
  double var3;
  double var4;
  int lo;
  int hi;
};

void *smoothing_thread(void *threadarg);
void *total_energy_thread(void *threadarg);
void *density_thread(void *threadarg);
void *pressure_thread(void *threadarg);
void *soundspeed_thread(void *threadarg);
void *CFL_thread(void *threadarg);
void *acceleration_thread(void *threadarg);
void *predictor_thread(void *threadarg);
void *corrector_thread(void *threadarg);

void create_smoothing_threads(struct universe *world, int iterations, int neighbours, double min_h,
			      double max_h, double *r, struct cell *tree, struct cell *root);
void create_density_threads(struct universe *world);
void create_pressure_threads(struct universe *world);
void create_soundspeed_threads(struct universe *world);
void create_CFL_threads(struct universe *world);
void create_acceleration_threads(struct universe *world);
void create_predictor_threads(struct universe *world);
void create_corrector_threads(struct universe *world);
void create_total_energy_threads(struct universe *world, double theta);

#endif
