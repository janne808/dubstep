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

__global__ void smoothing_length_iterator_kernel(dubfloat_t *r, dubfloat_t *origin, int *buffer, int *buffer_index, dubfloat_t *h, int n){
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
  //for(ii=lo;ii<hi;ii++){
  //}

  // wait for threads to finish
  __syncthreads();
}

void compute_smoothing_length_neighbours_cuda(struct universe *world, int iterations, int N_target){
  /* loop variables */
  int ii;
  int jj;
  int kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  dubfloat_t *r_in;
  dubfloat_t *h_in;

  /* device buffers */
  dubfloat_t *r_d;

  /* pointer to particle neighbour number vector */
  int *num_neighbours_in;

  /* maximum list length */
  int max_list_len;

  /* flattened neighbour list on device mem */
  int *flat_list_d;

  /* allocate particle displacement vector on device */
  cudaMalloc(&r_d, 3*world->num*sizeof(dubfloat_t));

  /* copy particle displacement vector on device */
  cudaMemcpy(r_d, world->r, 3*world->num*sizeof(dubfloat_t), cudaMemcpyHostToDevice);    

  /* flatten out neighbour list */
  /* search for max list size */
  max_list_len=0;
  for(ii=0;ii<world->num;ii++){
	if(world->neighbour_list[ii].max_size>max_list_len)
		max_list_len=world->neighbour_list[ii].max_size;
  }

  /* add slack to max size for new neighbours */
  max_list_len+=50;

  /* set up flat list to device memory */
  /* allocate list on device */
  cudaMalloc(&flat_list_d, world->num*max_list_len*sizeof(int));

  /* free device list */
  cudaFree(flat_list_d);

  /* free device particle displacement vector */
  cudaFree(r_d);
}

