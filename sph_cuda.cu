/* smoothed particle hydrodynamics cuda routines */

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
#include <math.h>
#include <cuda.h>

#include "dubstep.h"
#include "sph_cuda.h"

#define BLOCK_SIZE 8
#define NUM_BLOCKS 8

__global__ void smoothing_length_iterator_kernel(float *r, float *origin, int *buffer, int *buffer_index, float *h, int n){
  // loop index
  int ii;

  // thread index
  int idx=blockIdx.x*blockDim.x+threadIdx.x;

  // number of particles in thread slice
  int slice_len=n/(BLOCK_SIZE*NUM_BLOCKS);

  // low and high bound of thread slice
  int lo=slice_len*idx;
  int hi=slice_len*idx+slice_len;

  // last thread handles remaining particles
  if(idx==(BLOCK_SIZE*NUM_BLOCKS)){
    hi=n;
  }

  // find particle neighbours in every thread slice
  for(ii=lo;ii<hi;ii++){
  }

  // wait for threads to finish
  __syncthreads();
}

void compute_smoothing_length_neighbours_cuda(struct universe *world, int iterations, int N_target){
  /* loop variables */
  int ii,jj,kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *h_in;

  int *num_neighbours_in;

  /* host particle list buffer */
  int *buffer;

  /* device buffers */
  float *r_d;
  float *origin_d;

  int *buffer_d;
  int *buffer_index_d;

  float *h_d;
  float *h;

  int *buffer_index;

  float origin[3];

  // target for number of threads
  int num_threads=NUM_BLOCKS*BLOCK_SIZE;

  // compute execution configuration
  int blockSize=BLOCK_SIZE;
  int nBlocks=num_threads/blockSize;

  m=world->dim;  
  n=world->num;

  r_in=world->r2;
  h_in=world->h;

  num_neighbours_in=world->num_neighbours;

  // allocate particle list buffer on device
  cudaMalloc((void**)&buffer_d, num_threads*n*sizeof(int));
  cudaMalloc((void**)&buffer_index_d, num_threads*sizeof(int));

  // set up vector indeces
  buffer_index=(int*)malloc(num_threads*sizeof(int));
  for(ii=0;ii<num_threads;ii++)
    buffer_index[ii]=0;

  // allocate particle displacement vector on device
  cudaMalloc((void**)&r_d, m*n*sizeof(float));

  // allocate particle origin vector on device
  cudaMalloc((void**)&origin_d, m*sizeof(float));

  // copy particle displacement vector on device
  cudaMemcpy(r_d, r_in, m*n*sizeof(float), cudaMemcpyHostToDevice);

  // allocate smoothing length parameter on device
  cudaMalloc((void**)&h_d, sizeof(float));

  // allocate smoothing length parameter on host
  h=(float*)malloc(sizeof(float));

  // allocate particle list buffer
  buffer=(int*)malloc(n*sizeof(int));

  /* iterate towards optimum number of neighbours */
  for(ii=0;ii<n;ii++){
    // particle origin displacement vector
    origin[0]=(float)r_in[m*ii+0];
    origin[1]=(float)r_in[m*ii+1];
    origin[2]=(float)r_in[m*ii+2];

    *h=(float)h_in[ii];
    
    // copy smoothing length parameter to device
    cudaMemcpy(h_d, h, sizeof(float), cudaMemcpyHostToDevice);
      
    // copy particle origin vector on device
    cudaMemcpy(origin_d, &origin, m*sizeof(float), cudaMemcpyHostToDevice);

    // copy empty buffer vector index on device
    cudaMemcpy(buffer_index_d, buffer_index, num_threads*sizeof(int), cudaMemcpyHostToDevice);

    if(world->neighbour_list[ii].list){
      world->neighbour_list[ii].num=0;
      free(world->neighbour_list[ii].list);
    }

    // call kernel
    smoothing_length_iterator_kernel <<< nBlocks, blockSize >>> (r_d, origin_d, buffer_d, buffer_index_d, h_d, n);

    // copy smoothing length parameter to host
    cudaMemcpy(h, h_d, sizeof(float), cudaMemcpyDeviceToHost);

    h_in[ii]=(double)*h;

    num_neighbours_in[ii]=num_threads;
    
    world->neighbour_list[ii].num=num_threads;
    world->neighbour_list[ii].list=(int*)malloc(num_threads*sizeof(int));
    if(!world->neighbour_list[ii].list){
      printf("Out of memory: particle neighbour list not allocated.\n");
      exit(1);
    }
    memcpy(world->neighbour_list[ii].list, buffer, num_threads*sizeof(int));
  }

  // clean up
  free(buffer);
  free(h);
  cudaFree(h_d);
  cudaFree(r_d);
  cudaFree(origin_d);
  free(buffer_index);
  cudaFree(buffer_index_d);
  cudaFree(buffer_d);
}

