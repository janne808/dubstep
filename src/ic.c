/* initial condition generator routines */

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
#include "random.h"
#include "linear.h"
#include "ic.h"

void generate_glass(struct universe *world, dubfloat_t radius, dubfloat_t U_threshold, dubfloat_t epsilon, dubfloat_t G){
  /* loop variables */
  int ii;
  int jj;

  /* distance computation variables */
  dubfloat_t d;

  /* pointers to state vectors */
  dubfloat_t *r_in;
  dubfloat_t *m_in;

  /* total gravitational potential energy */
  dubfloat_t U;

  dubfloat_t rad;
  dubfloat_t theta;

  r_in=world->r;
  m_in=world->m;

  /* initial displacement */
  r_in[3*0+0]=0.0;
  r_in[3*0+1]=0.0;
  r_in[3*0+2]=0.0;

  for(ii=1;ii<world->num;ii++){
    //printf("Generating random displacement for particle %d...\n", ii);
    do{
      /* generate random displacement */
      rad=rejection_sampling();
      theta=2.0*PI*(dubfloat_t)rand()/RAND_MAX;
      r_in[3*ii+0]=rad*cos(theta);
      r_in[3*ii+1]=rad*sin(theta);
      r_in[3*ii+2]=50.0*boxmuller()*0.1;      

      //r_in[3*ii+0]=radius*boxmuller();
      //r_in[3*ii+1]=radius*boxmuller();
      //r_in[3*ii+2]=radius*boxmuller()*0.1;

      //do{
      //r_in[3*ii+0]=radius*(2.0*((dubfloat_t)rand()/RAND_MAX)-1.0);
      //r_in[3*ii+1]=radius*(2.0*((dubfloat_t)rand()/RAND_MAX)-1.0);
      //r_in[3*ii+2]=0.000001*radius*(2.0*((dubfloat_t)rand()/RAND_MAX)-1.0);
      //}while(sqrt(r_in[3*ii+0]*r_in[3*ii+0]+
      //	  r_in[3*ii+1]*r_in[3*ii+1]+
      //	  r_in[3*ii+2]*r_in[3*ii+2])>radius);

      /* calculate total gravitational energy for new displacement */
      /* unit is AU^3/(M_solar*yr^2)*M_solar^2/AU */
      U=0.0;
      for(jj=0;jj<ii;jj++){
	d=euclidean_distance(&r_in[ii*3],&r_in[jj*3],3);
	d=sqrt(d*d+epsilon*epsilon);
	U-=G*m_in[ii]*m_in[jj]*d;
      }
      //}while(U<U_threshold);
    }while(0);
  }
}

