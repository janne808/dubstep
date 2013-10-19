/* statistics code */

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
#include <string.h>
#include <math.h>

#include "dubstep.h"
#include "statistics.h"

#define CPUTIME_SMA_SIZE 16

struct sma_data *cputime_sma_init(){
  /* loop variables */
  int ii;

  /* data structs*/
  struct sma_data *data;

  /* allocate structure */
  data=(struct sma_data *)malloc(sizeof(struct sma_data));
  if(!data){
    printf("Out of memory: cputime SMA data struct not allocated.");
    exit(0);
  }

  /* allocate sample buffer */
  data->samples=(dubfloat_t *)malloc(CPUTIME_SMA_SIZE*sizeof(dubfloat_t));
  if(!data->samples){
    printf("Out of memory: cputime SMA sample buffer not allocated.");
    exit(0);
  }

  /* initialize structure variables */
  data->sum=0.0;

  for(ii=0;ii<CPUTIME_SMA_SIZE;ii++){
    data->samples[ii]=0.0;
  }

  data->write_point=0;

  /* return structure pointer */
  return data;
}

dubfloat_t cputime_sma_update(dubfloat_t input, struct sma_data *data){
  /* subtract last sample from the sum */
  data->sum-=data->samples[data->write_point]/(dubfloat_t)(CPUTIME_SMA_SIZE);

  /* insert input into sample buffer */
  data->samples[data->write_point++]=input;

  /* make sure write point is within buffer bounds */
  if(data->write_point>CPUTIME_SMA_SIZE-1)
    data->write_point-=CPUTIME_SMA_SIZE;

  /* add new sample to sum */
  data->sum+=input/(dubfloat_t)(CPUTIME_SMA_SIZE);

  /* return updated sum */
  return data->sum;
}

