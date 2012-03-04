/* dubstep */

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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dubstep.h"
#include "tree.h"
#include "sph.h"
#include "threads.h"

struct thread_data thread_data_array[NUM_THREADS+1];
struct thread_data2 thread_data_array2[NUM_THREADS+1];
struct thread_data3 thread_data_array3[NUM_THREADS+1];
struct thread_data4 thread_data_array4[NUM_THREADS+1];
struct thread_data5 thread_data_array5[NUM_THREADS+1];

void *smoothingThread(void *threadarg){
  struct thread_data4 *my_data;

  my_data=(struct thread_data4 *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_smoothing_length_neighbours(my_data->world, my_data->var1, my_data->var2, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *densityThread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_density(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *pressureThread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_pressure(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *soundspeedThread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_soundspeed(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *CFLThread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_cfl(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *accelerationThread(void *threadarg){
  struct thread_data2 *my_data;

  my_data=(struct thread_data2 *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_sph_acceleration(my_data->world, my_data->r, my_data->v, my_data->a, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *energyThread(void *threadarg){
  struct thread_data2 *my_data;

  my_data=(struct thread_data2 *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_internal_energy(my_data->world, my_data->a, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *integrationThread(void *threadarg){
  int nn;

  int m;
  int n;

  double r[3];
  double a[3];

  double *a_sph;
  double *a_tree;

  double *m_in;
  double *dt_CFL_in;

  double theta;
  double epsilon;
  double dt;

  struct cell *tree;
  struct universe *world;

  struct thread_data3 *my_data;

  my_data=(struct thread_data3 *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  world=my_data->world;

  m=world->dim;
  n=world->num;

  a_sph=world->a_sph;
  a_tree=world->a_tree;

  m_in=world->m;

  tree=world->tree;

  theta=my_data->var1;
  epsilon=my_data->var2;
  dt=world->sub_dt;
  dt_CFL_in=world->dt_CFL;

  /* integrate */
  for(nn=my_data->lo;nn<my_data->hi;nn++){
    /* displacement vector */
    r[0]=world->r[nn*m+0];
    r[1]=world->r[nn*m+1];
    r[2]=world->r[nn*m+2];

    world->u[nn]+=world->du[nn]*dt;

    /* approximate gravitational acceleration from tree */
    a[0]=0;
    a[1]=0;
    a[2]=0;

    forcerecurse(tree, &tree[0], r, a, world->G, theta, epsilon);

    a_tree[nn*m+0]=a[0];
    a_tree[nn*m+1]=a[1];
    a_tree[nn*m+2]=a[2];

    world->v[nn*m+0]+=a_tree[nn*m+0]*dt;
    world->v[nn*m+1]+=a_tree[nn*m+1]*dt;
    world->v[nn*m+2]+=a_tree[nn*m+2]*dt;

    world->v[nn*m+0]+=a_sph[nn*m+0]*dt;
    world->v[nn*m+1]+=a_sph[nn*m+1]*dt;
    world->v[nn*m+2]+=a_sph[nn*m+2]*dt;

    /* euler integration */
    world->r[nn*m+0]+=world->v[nn*m+0]*dt;
    world->r[nn*m+1]+=world->v[nn*m+1]*dt;
    world->r[nn*m+2]+=world->v[nn*m+2]*dt;
  }

  pthread_exit(NULL);
}

void *predictorThread(void *threadarg){
  int nn;

  int m;
  int n;

  double *a_sph;
  double *a_tree;

  double *m_in;
  double *dt_CFL_in;

  double theta;
  double epsilon;
  double dt;

  struct cell *tree;
  struct universe *world;

  struct thread_data5 *my_data;

  my_data=(struct thread_data5 *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  world=my_data->world;

  m=world->dim;
  n=world->num;

  a_sph=world->a_sph;
  a_tree=world->a_tree;

  m_in=world->m;

  tree=world->tree;

  theta=my_data->var1;
  epsilon=my_data->var2;
  dt=world->sub_dt;
  dt_CFL_in=world->dt_CFL;

  for(nn=my_data->lo;nn<my_data->hi;nn++){
    world->v2[nn*m+0]=world->v[nn*m+0]+world->a2[nn*m+0]*dt*0.5;
    world->v2[nn*m+1]=world->v[nn*m+1]+world->a2[nn*m+1]*dt*0.5;
    world->v2[nn*m+2]=world->v[nn*m+2]+world->a2[nn*m+2]*dt*0.5;

    world->r2[nn*m+0]=world->r[nn*m+0]+world->v2[nn*m+0]*dt*0.5;
    world->r2[nn*m+1]=world->r[nn*m+1]+world->v2[nn*m+1]*dt*0.5;
    world->r2[nn*m+2]=world->r[nn*m+2]+world->v2[nn*m+2]*dt*0.5;

    world->u2[nn]=world->u[nn]+world->du[nn]*dt*0.5;
    if(world->u2[nn]<1E-9)
      world->u2[nn]=1E-9;
  }

  pthread_exit(NULL);
}

void *correctorThread(void *threadarg){
  int nn;

  int m;
  int n;

  double r[3];
  double a[3];

  double *a_sph;
  double *a_tree;

  double *m_in;
  double *dt_CFL_in;

  double theta;
  double epsilon;
  double dt;

  struct cell *tree;
  struct universe *world;

  struct thread_data5 *my_data;

  my_data=(struct thread_data5 *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  world=my_data->world;

  m=world->dim;
  n=world->num;

  a_sph=world->a_sph;
  a_tree=world->a_tree;

  m_in=world->m;

  tree=world->tree;

  theta=my_data->var1;
  epsilon=my_data->var2;
  dt=world->sub_dt;
  dt_CFL_in=world->dt_CFL;

  /* integrate */
  for(nn=my_data->lo;nn<my_data->hi;nn++){
    /* displacement vector */
    r[0]=world->r2[nn*m+0];
    r[1]=world->r2[nn*m+1];
    r[2]=world->r2[nn*m+2];

    /* approximate gravitational acceleration from tree */
    a[0]=0;
    a[1]=0;
    a[2]=0;

    forcerecurse(tree, &tree[0], r, a, world->G, theta, epsilon);

    a_tree[nn*m+0]=a[0];
    a_tree[nn*m+1]=a[1];
    a_tree[nn*m+2]=a[2];

    world->a2[nn*m+0]=a_tree[nn*m+0]+a_sph[nn*m+0];
    world->a2[nn*m+1]=a_tree[nn*m+1]+a_sph[nn*m+1];
    world->a2[nn*m+2]=a_tree[nn*m+2]+a_sph[nn*m+2];

    world->v2[nn*m+0]=world->v[nn*m+0]+(world->a2[nn*m+0])*dt*0.5;
    world->v2[nn*m+1]=world->v[nn*m+1]+(world->a2[nn*m+1])*dt*0.5;
    world->v2[nn*m+2]=world->v[nn*m+2]+(world->a2[nn*m+2])*dt*0.5;

    world->r2[nn*m+0]=world->r[nn*m+0]+world->v2[nn*m+0]*dt*0.5;
    world->r2[nn*m+1]=world->r[nn*m+1]+world->v2[nn*m+1]*dt*0.5;
    world->r2[nn*m+2]=world->r[nn*m+2]+world->v2[nn*m+2]*dt*0.5;

    world->u2[nn]=world->u[nn]+world->du[nn]*dt*0.5;
    if(world->u2[nn]<1E-9)
      world->u2[nn]=1E-9;

    world->v[nn*m+0]=2*world->v2[nn*m+0]-world->v[nn*m+0];
    world->v[nn*m+1]=2*world->v2[nn*m+1]-world->v[nn*m+1];
    world->v[nn*m+2]=2*world->v2[nn*m+2]-world->v[nn*m+2];
    
    world->r[nn*m+0]=2*world->r2[nn*m+0]-world->r[nn*m+0];
    world->r[nn*m+1]=2*world->r2[nn*m+1]-world->r[nn*m+1];
    world->r[nn*m+2]=2*world->r2[nn*m+2]-world->r[nn*m+2];
    
    world->u[nn]=2*world->u2[nn]-world->u[nn];
    if(world->u[nn]<1E-9)
      world->u[nn]=1E-9;
  }

  pthread_exit(NULL);
}

void createSmoothingThreads(struct universe *world, int iterations, int neighbours){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array4[nn].thread_id=nn;
    thread_data_array4[nn].world=world;
    thread_data_array4[nn].var1=iterations;
    thread_data_array4[nn].var2=neighbours; 
    thread_data_array4[nn].lo=nn*thread_slice_num;
    thread_data_array4[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, smoothingThread, (void *) &thread_data_array4[nn]);

    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array4[nn].thread_id=nn;
    thread_data_array4[nn].world=world;
    thread_data_array4[nn].var1=iterations;
    thread_data_array4[nn].var2=neighbours;
    thread_data_array4[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array4[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, densityThread, (void *) &thread_data_array4[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }
}

void createDensityThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=1.0;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, densityThread, (void *) &thread_data_array[nn]);

    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }

    num_join_threads++;
  }

  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=1.0;
    thread_data_array[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, densityThread, (void *) &thread_data_array[nn]);

    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }
}

void createPressureThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=world->gamma;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, pressureThread, (void *) &thread_data_array[nn]);

    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }

  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=world->gamma;
    thread_data_array[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, pressureThread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }

    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }
}

void createSoundspeedThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=world->gamma;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, soundspeedThread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=world->gamma;
    thread_data_array[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, soundspeedThread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }
}

void createCFLThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=0.3;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, CFLThread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=0.3;
    thread_data_array[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, CFLThread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }
}

void createEnergyThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array2[nn].thread_id=nn;
    thread_data_array2[nn].world=world;
    thread_data_array2[nn].a=world->a_sph;
    thread_data_array2[nn].lo=nn*thread_slice_num;
    thread_data_array2[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, energyThread, (void *) &thread_data_array2[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array2[nn].thread_id=nn;
    thread_data_array2[nn].world=world;
    thread_data_array2[nn].a=world->a_sph;
    thread_data_array2[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array2[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, energyThread, (void *) &thread_data_array2[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }
}

void createAccelerationThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array2[nn].thread_id=nn;
    thread_data_array2[nn].world=world;
    thread_data_array2[nn].r=world->r2;
    thread_data_array2[nn].v=world->v2;
    thread_data_array2[nn].a=world->a_sph;
    thread_data_array2[nn].lo=nn*thread_slice_num;
    thread_data_array2[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, accelerationThread, (void *) &thread_data_array2[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array2[nn].thread_id=nn;
    thread_data_array2[nn].world=world;
    thread_data_array2[nn].r=world->r2;
    thread_data_array2[nn].v=world->v2;
    thread_data_array2[nn].a=world->a_sph;
    thread_data_array2[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array2[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, soundspeedThread, (void *) &thread_data_array2[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }
}

void createIntegrationThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array3[nn].thread_id=nn;
    thread_data_array3[nn].world=world;
    thread_data_array3[nn].a_sph=world->a_sph;
    thread_data_array3[nn].var1=world->epsilon;
    thread_data_array3[nn].var2=world->theta;
    thread_data_array3[nn].var3=world->dt;
    thread_data_array3[nn].lo=nn*thread_slice_num;
    thread_data_array3[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, integrationThread, (void *) &thread_data_array3[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array3[nn].thread_id=nn;
    thread_data_array3[nn].world=world;
    thread_data_array3[nn].a_sph=world->a_sph;
    thread_data_array3[nn].var1=world->epsilon;
    thread_data_array3[nn].var2=world->theta;
    thread_data_array3[nn].var3=world->dt;
    thread_data_array3[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array3[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, integrationThread, (void *) &thread_data_array3[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }      
}

void createPredictorThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array5[nn].thread_id=nn;
    thread_data_array5[nn].world=world;
    thread_data_array5[nn].a_sph=world->a_sph;
    thread_data_array5[nn].var1=world->epsilon;
    thread_data_array5[nn].var2=world->theta;
    thread_data_array5[nn].var3=world->dt;
    thread_data_array5[nn].lo=nn*thread_slice_num;
    thread_data_array5[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, predictorThread, (void *) &thread_data_array5[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array5[nn].thread_id=nn;
    thread_data_array5[nn].world=world;
    thread_data_array5[nn].a_sph=world->a_sph;
    thread_data_array5[nn].var1=world->epsilon;
    thread_data_array5[nn].var2=world->theta;
    thread_data_array5[nn].var3=world->dt;
    thread_data_array5[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array5[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, predictorThread, (void *) &thread_data_array5[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }      
}

void createCorrectorThreads(struct universe *world){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array5[nn].thread_id=nn;
    thread_data_array5[nn].world=world;
    thread_data_array5[nn].a_sph=world->a_sph;
    thread_data_array5[nn].var1=world->epsilon;
    thread_data_array5[nn].var2=world->theta;
    thread_data_array5[nn].var3=world->dt;
    thread_data_array5[nn].lo=nn*thread_slice_num;
    thread_data_array5[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, correctorThread, (void *) &thread_data_array5[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                world->num);
    thread_data_array5[nn].thread_id=nn;
    thread_data_array5[nn].world=world;
    thread_data_array5[nn].a_sph=world->a_sph;
    thread_data_array5[nn].var1=world->epsilon;
    thread_data_array5[nn].var2=world->theta;
    thread_data_array5[nn].var3=world->dt;
    thread_data_array5[nn].lo=(nn-1)*thread_slice_num;
    thread_data_array5[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, correctorThread, (void *) &thread_data_array5[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  /* join threads */
  pthread_attr_destroy(&attr);
  for(nn=0;nn<num_join_threads;nn++){
    thread_rc=pthread_join(threads[nn], &thread_status);
    
    if(thread_rc){
      printf("ERROR: pthread_join() returned %d.\n", thread_rc);
      exit(-1);
    }
  }      
}
