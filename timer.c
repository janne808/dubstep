/* timer routines */

/*
 *  (C) 2012 Janne Heikkarainen <janne.heikkarainen@tut.fi>
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

#include <time.h>

#include "timer.h"

/* compute time difference in seconds and nanoseconds */
void timediff(struct timespec start, struct timespec end, struct timespec *out){
  /* compute time difference */
  /* handle nsec overflow as needed */
  if((end.tv_nsec-start.tv_nsec)<0){
    out->tv_nsec=1000000000-start.tv_nsec+end.tv_nsec;
    out->tv_sec=end.tv_sec-start.tv_sec-1;
  }
  else{
    out->tv_nsec=end.tv_nsec-start.tv_nsec;
    out->tv_sec=end.tv_sec-start.tv_sec;
  }
}
