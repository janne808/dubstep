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

void *smoothingThread(void *threadarg);
void *densityThread(void *threadarg);
void *pressureThread(void *threadarg);
void *soundspeedThread(void *threadarg);
void *CFLThread(void *threadarg);
void *energyThread(void *threadarg);
void *accelerationThread(void *threadarg);
void *integrationThread(void *threadarg);
void *predictorThread(void *threadarg);
void *correctorThread(void *threadarg);

void createSmoothingThreads(struct universe *world, int iterations, int neighbours);
void createDensityThreads(struct universe *world);
void createPressureThreads(struct universe *world);
void createSoundspeedThreads(struct universe *world);
void createCFLThreads(struct universe *world);
void createEnergyThreads(struct universe *world);
void createAccelerationThreads(struct universe *world);
void createIntegrationThreads(struct universe *world);
void createPredictorThreads(struct universe *world);
void createCorrectorThreads(struct universe *world);

#endif
