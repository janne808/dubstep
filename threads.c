/* thread routines */

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
struct thread_data6 thread_data_array6[NUM_THREADS+1];

void *smoothing_thread(void *threadarg){
  struct thread_data6 *my_data;

  my_data=(struct thread_data6 *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

#if (defined ADAPTIVE_SMOOTHING)&&ADAPTIVE_SMOOTHING
  compute_smoothing_length_tree(my_data->world, my_data->var1, my_data->var2, my_data->var3, my_data->var4,
  				my_data->r, my_data->tree, my_data->root, my_data->lo, my_data->hi);
#else
  compute_constant_smoothing_length_tree(my_data->world, my_data->var1, my_data->var2, my_data->var3, my_data->var4,
  					 my_data->r, my_data->tree, my_data->root, my_data->lo, my_data->hi);
#endif

  pthread_exit(NULL);
}

void *total_energy_thread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_total_energy(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *density_thread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_density(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *pressure_thread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_pressure(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *soundspeed_thread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_soundspeed(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *CFL_thread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_cfl(my_data->world, my_data->var, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *timebin_thread(void *threadarg){
  struct thread_data *my_data;

  my_data=(struct thread_data *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  update_time_bins(my_data->world, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *acceleration_thread(void *threadarg){
  struct thread_data2 *my_data;

  my_data=(struct thread_data2 *) threadarg;
  //printf("Executing thread %d, slice %d to %d.\n", my_data->thread_id, my_data->lo, my_data->hi);

  compute_internal_energy_and_acceleration(my_data->world, my_data->r, my_data->v, my_data->a, my_data->lo, my_data->hi);

  pthread_exit(NULL);
}

void *predictor_thread(void *threadarg){
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

void *corrector_thread(void *threadarg){
  int nn;

  int m;
  int n;

  double r[3];
  double a[3];

  double *a_sph;
  double *a_tree;

  double *m_in;
  double *h_in;

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
  h_in=world->h;

  tree=world->tree;

  theta=my_data->var1;
  epsilon=my_data->var2;
  dt_CFL_in=world->dt_CFL;

  /* integrate */
  for(nn=my_data->lo;nn<my_data->hi;nn++){
    /* displacement vector */
    r[0]=world->r2[nn*m+0];
    r[1]=world->r2[nn*m+1];
    r[2]=world->r2[nn*m+2];

    /* kick particles with individual time steps */
    if(world->kick[nn]){
      /* compute individual time step for kick */
      dt=world->dt/pow(2,world->time_bin[nn]);

      /* approximate gravitational acceleration from tree */
      a[0]=0;
      a[1]=0;
      a[2]=0;

#if (defined ADAPTIVE_SOFTENING)&&ADAPTIVE_SOFTENING
      force_walk(world, tree, &tree[0], r, a, world->G, theta, h_in[nn]);
#else
      //force_walk(world, tree, &tree[0], r, a, world->G, theta, world->epsilon);
      direct_summation(world, r, a, world->G);
#endif

      a_tree[nn*m+0]=a[0];
      a_tree[nn*m+1]=a[1];
      a_tree[nn*m+2]=a[2];

      world->a2[nn*m+0]=a_tree[nn*m+0]+a_sph[nn*m+0];
      world->a2[nn*m+1]=a_tree[nn*m+1]+a_sph[nn*m+1];
      world->a2[nn*m+2]=a_tree[nn*m+2]+a_sph[nn*m+2];

      world->v2[nn*m+0]=world->v[nn*m+0]+(world->a2[nn*m+0])*dt*0.5;
      world->v2[nn*m+1]=world->v[nn*m+1]+(world->a2[nn*m+1])*dt*0.5;
      world->v2[nn*m+2]=world->v[nn*m+2]+(world->a2[nn*m+2])*dt*0.5;

      world->v[nn*m+0]=2*world->v2[nn*m+0]-world->v[nn*m+0];
      world->v[nn*m+1]=2*world->v2[nn*m+1]-world->v[nn*m+1];
      world->v[nn*m+2]=2*world->v2[nn*m+2]-world->v[nn*m+2];    
    }
    else{
      world->v2[nn*m+0]=world->v[nn*m+0];
      world->v2[nn*m+1]=world->v[nn*m+1];
      world->v2[nn*m+2]=world->v[nn*m+2];
    }
    
    /* integrate position and internal energy with smallest time step */
    dt=world->sub_dt;

    world->r2[nn*m+0]=world->r[nn*m+0]+world->v2[nn*m+0]*dt*0.5;
    world->r2[nn*m+1]=world->r[nn*m+1]+world->v2[nn*m+1]*dt*0.5;
    world->r2[nn*m+2]=world->r[nn*m+2]+world->v2[nn*m+2]*dt*0.5;

    world->u2[nn]=world->u[nn]+world->du[nn]*dt*0.5;
    if(world->u2[nn]<1E-9)
      world->u2[nn]=1E-9;

    world->r[nn*m+0]=2*world->r2[nn*m+0]-world->r[nn*m+0];
    world->r[nn*m+1]=2*world->r2[nn*m+1]-world->r[nn*m+1];
    world->r[nn*m+2]=2*world->r2[nn*m+2]-world->r[nn*m+2];
    
    world->u[nn]=2*world->u2[nn]-world->u[nn];
    if(world->u[nn]<1E-9)
      world->u[nn]=1E-9;
  }

  pthread_exit(NULL);
}

void create_smoothing_threads(struct universe *world, int iterations, int neighbours, double min_h,
			      double max_h, double *r, struct cell *tree, struct cell *root){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array6[nn].thread_id=nn;
    thread_data_array6[nn].world=world;
    thread_data_array6[nn].var1=min_h;
    thread_data_array6[nn].var2=max_h;
    thread_data_array6[nn].var3=iterations;
    thread_data_array6[nn].var4=neighbours; 
    thread_data_array6[nn].r=r;
    thread_data_array6[nn].tree=tree;
    thread_data_array6[nn].root=root;
    thread_data_array6[nn].lo=nn*thread_slice_num;
    thread_data_array6[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, smoothing_thread, (void *) &thread_data_array6[nn]);

    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array6[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array6[nn-1].hi,
    //                                                world->num);
    thread_data_array6[nn].thread_id=nn;
    thread_data_array6[nn].world=world;
    thread_data_array6[nn].var1=min_h;
    thread_data_array6[nn].var2=max_h;
    thread_data_array6[nn].var3=iterations;
    thread_data_array6[nn].var4=neighbours; 
    thread_data_array6[nn].r=r;
    thread_data_array6[nn].tree=tree;
    thread_data_array6[nn].root=root;
    thread_data_array6[nn].lo=thread_data_array6[nn-1].hi;
    thread_data_array6[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, smoothing_thread, (void *) &thread_data_array6[nn]);
    
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

void create_density_threads(struct universe *world){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=1.0;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, density_thread, (void *) &thread_data_array[nn]);

    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }

    num_join_threads++;
  }

  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array[nn-1].hi,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=1.0;
    thread_data_array[nn].lo=thread_data_array[nn-1].hi;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, density_thread, (void *) &thread_data_array[nn]);

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

void create_pressure_threads(struct universe *world){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=world->gamma;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, pressure_thread, (void *) &thread_data_array[nn]);

    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }

  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array[nn-1].hi,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=world->gamma;
    thread_data_array[nn].lo=thread_data_array[nn-1].hi;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, pressure_thread, (void *) &thread_data_array[nn]);
    
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

void create_soundspeed_threads(struct universe *world){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=world->gamma;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, soundspeed_thread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array[nn-1].hi,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=world->gamma;
    thread_data_array[nn].lo=thread_data_array[nn-1].hi;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, soundspeed_thread, (void *) &thread_data_array[nn]);
    
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

void create_CFL_threads(struct universe *world){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=0.3;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, CFL_thread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array[nn-1].hi,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=0.3;
    thread_data_array[nn].lo=thread_data_array[nn-1].hi;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, CFL_thread, (void *) &thread_data_array[nn]);
    
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

void create_timebin_threads(struct universe *world){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=0.3;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, timebin_thread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array[nn-1].hi,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=0.3;
    thread_data_array[nn].lo=thread_data_array[nn-1].hi;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, timebin_thread, (void *) &thread_data_array[nn]);
    
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

void create_acceleration_threads(struct universe *world){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array2[nn].thread_id=nn;
    thread_data_array2[nn].world=world;
    thread_data_array2[nn].r=world->r2;
    thread_data_array2[nn].v=world->v2;
    thread_data_array2[nn].a=world->a_sph;
    thread_data_array2[nn].lo=nn*thread_slice_num;
    thread_data_array2[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, acceleration_thread, (void *) &thread_data_array2[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array2[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array2[nn-1].hi,
    //                                                world->num);
    thread_data_array2[nn].thread_id=nn;
    thread_data_array2[nn].world=world;
    thread_data_array2[nn].r=world->r2;
    thread_data_array2[nn].v=world->v2;
    thread_data_array2[nn].a=world->a_sph;
    thread_data_array2[nn].lo=thread_data_array2[nn-1].hi;
    thread_data_array2[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, acceleration_thread, (void *) &thread_data_array2[nn]);
    
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

void create_predictor_threads(struct universe *world){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
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
    thread_rc=pthread_create(&threads[nn], &attr, predictor_thread, (void *) &thread_data_array5[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array5[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array5[nn-1].hi,
    //                                                world->num);
    thread_data_array5[nn].thread_id=nn;
    thread_data_array5[nn].world=world;
    thread_data_array5[nn].a_sph=world->a_sph;
    thread_data_array5[nn].var1=world->epsilon;
    thread_data_array5[nn].var2=world->theta;
    thread_data_array5[nn].var3=world->dt;
    thread_data_array5[nn].lo=thread_data_array5[nn-1].hi;
    thread_data_array5[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, predictor_thread, (void *) &thread_data_array5[nn]);
    
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

void create_corrector_threads(struct universe *world){
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
  for(nn=0;nn<NUM_THREADS-1;nn++){
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
    thread_rc=pthread_create(&threads[nn], &attr, corrector_thread, (void *) &thread_data_array5[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array5[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array5[nn-1].hi,
    //                                                world->num);
    thread_data_array5[nn].thread_id=nn;
    thread_data_array5[nn].world=world;
    thread_data_array5[nn].a_sph=world->a_sph;
    thread_data_array5[nn].var1=world->epsilon;
    thread_data_array5[nn].var2=world->theta;
    thread_data_array5[nn].var3=world->dt;
    thread_data_array5[nn].lo=thread_data_array5[nn-1].hi;
    thread_data_array5[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, corrector_thread, (void *) &thread_data_array5[nn]);
    
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

void create_total_energy_threads(struct universe *world, double theta){
  /* posix thread variables */
  int thread_slice_num;
  int num_join_threads;
  int thread_rc;
  pthread_t threads[NUM_THREADS+1];
  pthread_attr_t attr;
  void *thread_status;
  
  /* loop variables */
  int nn;

  /* init gravitational potential and kinetic energy */
  world->u_grav=0;
  world->u_kin=0;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_slice_num=world->num/NUM_THREADS;
  num_join_threads=0;
  for(nn=0;nn<NUM_THREADS-1;nn++){
    //printf("Creating thread %d, slice %d to %d.\n", nn, nn*thread_slice_num,
    //                                                nn*thread_slice_num+thread_slice_num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=theta;
    thread_data_array[nn].lo=nn*thread_slice_num;
    thread_data_array[nn].hi=nn*thread_slice_num+thread_slice_num;
    thread_rc=pthread_create(&threads[nn], &attr, total_energy_thread, (void *) &thread_data_array[nn]);
    
    if(thread_rc){
      printf("ERROR: pthread_create() returned %d.\n", thread_rc);
      exit(-1);
    }
    
    num_join_threads++;
  }
  
  if(thread_data_array[nn-1].hi<world->num){
    //printf("Creating thread %d, slice %d to %d.\n", nn, thread_data_array[nn-1].hi,
    //                                                world->num);
    thread_data_array[nn].thread_id=nn;
    thread_data_array[nn].world=world;
    thread_data_array[nn].var=theta;
    thread_data_array[nn].lo=thread_data_array[nn-1].hi;
    thread_data_array[nn].hi=world->num;
    thread_rc=pthread_create(&threads[nn], &attr, total_energy_thread, (void *) &thread_data_array[nn]);
    
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

