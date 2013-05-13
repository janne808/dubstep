/* dubstep routines */

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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

#if ENABLE_GUI
#include <SDL.h>
#include <SDL_opengl.h>
#include <tiffio.h>
#endif

#include "dubstep.h"
#include "tree.h"
#include "sph.h"
#include "threads.h"
#include "timer.h"
#include "random.h"
#include "ic.h"
#include "linear.h"

#if ENABLE_GUI
void colormap(double val, struct color *col){
  int nn;

  double x_map[6]={0.0, 0.2, 0.45, 0.7, 0.85, 1.0};
  double r_map[6]={0.0, 0.0, 0.0, 1.0, 1.0, 0.65};
  double g_map[6]={0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
  double b_map[6]={0.65, 1.0, 1.0, 0.0, 0.0, 0.0};

  /* crop value to fit the map */
  if(val>0.99)
    val=0.99;

  /* linearly interpolate the value from the colormap */
  for(nn=0;nn<(6-1);nn++){
    if(val>x_map[nn]&&val<x_map[nn+1]){
      col->r=r_map[nn]+(val-x_map[nn])*(r_map[nn+1]-r_map[nn])/(x_map[nn+1]-x_map[nn]);
      col->g=g_map[nn]+(val-x_map[nn])*(g_map[nn+1]-g_map[nn])/(x_map[nn+1]-x_map[nn]);
      col->b=b_map[nn]+(val-x_map[nn])*(b_map[nn+1]-b_map[nn])/(x_map[nn+1]-x_map[nn]);
      break;
    }
  }
}

void writeframe(char* path){
  TIFF *file;
  GLubyte *image, *p;
  int ii;

  int width=640, height=480;

  file=TIFFOpen(path,"w");

  if(file){
    image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);

    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
    TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
    TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
    TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(file, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
    TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, "");

    p = image;
    for (ii = height - 1; ii >= 0; ii--) {
      if (TIFFWriteScanline(file, p, ii, 0) < 0) {
	free(image);
	TIFFClose(file);
	printf("Error writing TIFF file.\n");
	exit(1);
      }
      p += width * sizeof(GLubyte) * 3;
    }
    free(image);
    TIFFClose(file);
  }
}

#endif

int main(int argc, char *argv[])
{
  /* loop indices */
  int ii=0, nn=0, tt=0;
    
  /* matrix dimension variables */
  int m, n;
    
  /* tree cell array */
  struct cell *tree=0;

  /* physical variables */
  double *r, *a_tree, *a_sph;

  /* computational model parameters */

  /* gravitation constant normalized to astronomical units */
  /* i.e. mass is normalized to solar mass, time to 1 year and */
  /* distance to AU */
  double G=4.0*PI*PI;

  /* tree cell opening length parameter */
  double theta=0.85;

  /* maximum timestep */
  double dt=0.085;

  /* plummer gravitational softening factor*/
  double epsilon=SOFTENING_FACTOR;

  /* initial smoothing length */
  double h=MAX_SMOOTH_LEN;

  /* artificial viscosity parameters */
  double alpha=1.0;
  double beta=1.0;

  /* adiabatic exponent*/
  double gamma=1.4;

  /* time bin handling variables */
  int max_bin;

#if ENABLE_GUI
  /* SDL variables */
  const SDL_VideoInfo *info;
  int bpp;
  int screen_width=640, screen_height=480;

  SDL_Event *event;
#endif

  /* run flag for main loop */
  int run;
  int calc;

  /* time measurement variables */
  struct timespec treetime;
  struct timespec int_time;
  struct timespec int2_time;
  struct timespec sph_time;

  struct timespec time1, time2;

  /* particle neighbour statistics */
  int max_N;
  int min_N;
  double avg_N;

#if (defined ENERGY_PROFILING)&&ENERGY_PROFILING
  /* variables for energy conservation profiling */
  double total_u;
  double total_u2;
#endif

#if ENABLE_GUI
  /* UI variables */
  double rot[3];
  double rot2[3];
  int mouse_x,mouse_y;
  int mouse_dx,mouse_dy;
  int buttonstate;
  double zoom=1/5.0;
  double zoom2;

  double density;
  double max_density;

  double pressure;
  double max_pressure;

  struct color *col=0;

  int tiff_frame=0;
#endif

#if ENABLE_GUI
  char filename[128];
#endif

  /* world structure */
  struct universe *world=0;

#if ENABLE_GUI
  event=(SDL_Event*)malloc(sizeof(SDL_Event));
  if(!event){
    printf("Out of memory: SDL event not allocated.\n");
    exit(1);
  }      

  /* initialize graphics */
  if(SDL_Init(SDL_INIT_VIDEO)<0){
    printf("Unable to initialize SDL.\n");
    SDL_Quit();
    exit(1);
  }

  info=SDL_GetVideoInfo();
  if(!info){
    printf("Video query failed.\n");
    SDL_Quit();
    exit(1);
  }

  bpp=info->vfmt->BitsPerPixel;

  // set bits for red: (5 = 5bits for red channel)
  SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5);
  // set bits for green:
  SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 5);
  // set bits for blue:
  SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5);
  // colour depth:
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
  // you want it double buffered?
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, true);

  if(SDL_SetVideoMode(screen_width, screen_height, bpp,
		      SDL_OPENGL | SDL_SWSURFACE) == 0){
    printf("Unable to set video mode.\n");
    SDL_Quit();
    exit(1);
  }

  /* set up opengl */
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glViewport(0, 0, screen_width, screen_height);
  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-10.0f, 10.0f, 10.0f, -10.0f, -10.0f, 10.0f);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  col=(struct color*)malloc(sizeof(struct color));
  if(!col){
    printf("Out of memory: col not allocated.\n");
    exit(1);
  }
#endif

  /* allocate tree and vectors */
  tree=(struct cell*)malloc(100000*sizeof(struct cell));
  if(!tree){
    printf("Out of memory: tree not allocated.\n");
    exit(1);
  }      
      
  /* allocate world */
  world=(struct universe*)malloc(sizeof(struct universe));
  if(!world){
    printf("Out of memory: world not allocated.\n");
    exit(1);
  }

  /* set up world structure */
  world->tree=tree;

  world->G=G;

  world->dim=3;
  world->num=NUM;

  m=world->dim;
  n=world->num;

  world->dt=dt;
  world->theta=theta;
  world->epsilon=epsilon;
  world->gamma=gamma;
  world->alpha=alpha;
  world->beta=beta;

  world->neighbour_list=(struct particlelist*)malloc(n*sizeof(struct particlelist));
  /* allocate neighbour list for initial neighbours */
  /* neighbourhood searcher reallocates more memory when necessary */
  for(ii=0;ii<world->num;ii++){
    world->neighbour_list[ii].list=(int*)malloc(100*sizeof(int));
    if(!world->neighbour_list[ii].list){
      printf("Out of memory: particle neighbour list not allocated.\n");
      exit(1);
    }
    world->neighbour_list[ii].max_size=100;
  }

  world->m=(double*)malloc(n*sizeof(double));

  world->r=(double*)malloc(m*n*sizeof(double));
  world->r2=(double*)malloc(m*n*sizeof(double));
  world->v=(double*)malloc(m*n*sizeof(double));
  world->v2=(double*)malloc(m*n*sizeof(double));
  world->a=(double*)malloc(m*n*sizeof(double));
  world->a2=(double*)malloc(m*n*sizeof(double));

  world->rho=(double*)malloc(n*sizeof(double));
  world->p=(double*)malloc(n*sizeof(double));
  world->c=(double*)malloc(n*sizeof(double));

  world->u=(double*)malloc(n*sizeof(double));
  world->u2=(double*)malloc(n*sizeof(double));
  world->du=(double*)malloc(n*sizeof(double));
  world->du2=(double*)malloc(n*sizeof(double));

  world->del_rho=(double*)malloc(n*sizeof(double));

  world->h=(double*)malloc(n*sizeof(double));

  world->dt_CFL=(double*)malloc(n*sizeof(double));
  world->kick=(int*)malloc(n*sizeof(int));
  world->kick_list=(int*)malloc(n*sizeof(int));
  world->time_bin=(int*)malloc(n*sizeof(int));

  world->cellindex=(int*)malloc(n*sizeof(int));

  r=(double*)malloc(m*sizeof(double));
  a_tree=(double*)malloc(m*n*sizeof(double));
  a_sph=(double*)malloc(m*n*sizeof(double));

  world->num_neighbours=(int*)malloc(n*sizeof(int));

  world->timediv=(int)(pow(2,3));
  world->sub_dt=world->dt;

  /* init world time */
  world->time=0.0;

  /* init pseudorandom number generator */
  srand(32423423);

  /* initial mass and velocities */
  for(ii=0;ii<n;ii++){
    world->neighbour_list[ii].num=0;
    world->num_neighbours[ii]=0;
    world->cellindex[ii]=0;
    world->h[ii]=h;
    world->dt_CFL[ii]=world->sub_dt;
    world->kick[ii]=1;
    world->time_bin[ii]=0;
    world->m[ii]=0.25/(double)(n);
    world->v[ii*m+0]=0;
    world->v[ii*m+1]=0;
    world->v[ii*m+2]=0;
    world->v2[ii*m+0]=world->v[ii*m+0];
    world->v2[ii*m+1]=world->v[ii*m+1];
    world->v2[ii*m+2]=world->v[ii*m+2];
    world->a[ii*m+0]=0;
    world->a[ii*m+1]=0;
    world->a[ii*m+2]=0;
    world->a2[ii*m+0]=world->a[ii*m+0];
    world->a2[ii*m+1]=world->a[ii*m+1];
    world->a2[ii*m+2]=world->a[ii*m+2];
    a_sph[ii*m+0]=0;
    a_sph[ii*m+1]=0;
    a_sph[ii*m+2]=0;
    a_tree[ii*m+0]=0;
    a_tree[ii*m+1]=0;
    a_tree[ii*m+2]=0;
  }

  world->a_sph=a_sph;
  world->a_tree=a_tree;

  /* initial thermal energy */
  for(ii=0;ii<n;ii++){
    world->u[ii]=0.2;
    world->u2[ii]=world->u[ii];
  }

  /* initial displacement */
  generate_glass(world, 15.0, -0.01*G, epsilon, G);

  for(ii=0;ii<n;ii++){
    double rr;
    double vv;

    world->r2[ii*3+0]=world->r[ii*3+0];
    world->r2[ii*3+1]=world->r[ii*3+1];
    world->r2[ii*3+2]=world->r[ii*3+2];

    rr=euclidean_norm(&world->r[ii*3],3)+1.0E-8;

    vv=sqrt(2.0*G*0.05125*2.0/rr/2.0);
    world->v[ii*3+0]=-vv*world->r[ii*3+1]/rr;
    world->v[ii*3+1]=vv*world->r[ii*3+0]/rr;
    world->v[ii*3+2]=0;    

    world->v2[ii*3+0]=world->v[ii*3+0];
    world->v2[ii*3+1]=world->v[ii*3+1];
    world->v2[ii*3+2]=world->v[ii*3+2];    
  }

  /* compute initial sph variables */

  /* create threads for smoothing length interation and
     interacting particle list generation */
  //create_smoothing_threads(world, 10, 25);
      
  /* initialize tree root parameters */
  init_treeroot(tree, world, world->r2);
    
  /* form tree */
  branch_recurse(world, tree, &tree[0], world->cellindex);
      
  /* init kick list for SPH computation */
  init_kick_list(world);

  /* create threads for parallel tree smoothing length iterators */
  create_smoothing_threads(world, 1, 25, MIN_SMOOTH_LEN, MAX_SMOOTH_LEN, world->r2, tree, &tree[0]);

  /* serial tree smoothing length iterator */
  //compute_smoothing_length_tree(world, MIN_SMOOTH_LEN, MAX_SMOOTH_LEN, 10, 25, world->r2, tree, &tree[0], 0, world->num);

  /* create threads for density computation */
  create_density_threads(world);

  /* create threads for pressure computation */
  create_pressure_threads(world);

  /* create threads for sound speed computation */
  create_soundspeed_threads(world);

  /* filter initial velocity field */
  smooth_velocity_field(world, 0, world->num);

  /* create threads for CFL computation */
  create_CFL_threads(world);

  /* create threads for hydrodynamic acceleration and internal energy computation */
  create_acceleration_threads(world);

  /* filter initial energy field */
  smooth_energy_field(world, 0, world->num);

  /* init particle time bins */
  for(nn=0;nn<world->num;nn++){
    for(ii=0;ii<world->num;ii++){
      if(world->dt_CFL[nn]>world->dt/pow(2,ii)){
	world->time_bin[nn]=ii;
	break;
      }
    }
    //printf("%d\t", world->time_bin[nn]);
  }
  //printf("\n");

  /* determine minimum timestep to take from maximum bin */
  max_bin=0;
  for(nn=0;nn<world->num;nn++){
    if(world->time_bin[nn]>max_bin)
      max_bin=world->time_bin[nn];
  }
  world->sub_dt=world->dt/pow(2,max_bin+1);

  /* integrate corrector */
  create_corrector_threads(world);

  /* initialize run flag and run the main loop */
  tt=0;
  run=1;
  calc=1;

#if ENABLE_GUI
  rot[0]=0;
  rot[1]=0;
  rot2[0]=rot[0];
  rot2[1]=rot[1];
  zoom2=0;
#endif

  while(run){
#if ENABLE_GUI
    /* check and handle events */
    while(SDL_PollEvent(event)){
      switch(event->type){
      case SDL_MOUSEMOTION:
	if(calc){
	  mouse_x=event->motion.x;
	  mouse_y=event->motion.y;
	  mouse_dx=mouse_x;
	  mouse_dy=mouse_y;
	}
	else{ 
	  if(buttonstate==SDL_BUTTON_LEFT){
	    mouse_x=event->motion.x-mouse_dx;
	    mouse_y=event->motion.y-mouse_dy;
	    rot2[0]=0.3*mouse_x;
	    rot2[1]=0.3*mouse_y;
	  }else if(buttonstate==SDL_BUTTON_RIGHT){
	    mouse_y=event->motion.y-mouse_dy;
	    zoom*=(1+0.0001*mouse_y);
	  }
	}
	break;
      case SDL_MOUSEBUTTONDOWN:
	if(calc){
	  mouse_dx=mouse_x;
	  mouse_dy=mouse_y;
	  calc=0;
	  buttonstate=(int)(event->button.button);
	}
	break;
      case SDL_MOUSEBUTTONUP:
	calc=1;
	if(buttonstate==SDL_BUTTON_LEFT){
	  rot[0]+=rot2[0];
	  rot[1]+=rot2[1];
	  rot2[0]=0;
	  rot2[1]=0;
	}else if(buttonstate==SDL_BUTTON_RIGHT){
	  //zoom+=zoom2;
	  //zoom2=0;
	}
	buttonstate=0;
	break;
      case SDL_QUIT:
	run=0;
	break;
      default:
	break;
      }
    }
#endif

    if(calc){
      
      /* timer start */
      clock_gettime(CLOCK_MONOTONIC, &time1);

      /* integrate predictor */
      create_predictor_threads(world);

      /* timer stop */
      clock_gettime(CLOCK_MONOTONIC, &time2);

      /* compute time */
      timediff(time1, time2, &int_time);

      /* compute sph variables */

      /* create threads for smoothing length interation and
	 interacting particle list generation */

      /* timer start */
      clock_gettime(CLOCK_MONOTONIC, &time1);

      /* initialize tree root parameters */
      init_treeroot(tree, world, world->r2);
    
      /* form tree */
      branch_recurse(world, tree, &tree[0], world->cellindex);

      /* timer stop */
      clock_gettime(CLOCK_MONOTONIC, &time2);

      /* compute time */
      timediff(time1, time2, &treetime);

      /* timer start */
      clock_gettime(CLOCK_MONOTONIC, &time1);

      /* update kick list for SPH computation */
      update_kick_list(world);

      /* serial tree smoothing length iterator */
      //compute_smoothing_length_tree(world, h, 1, 25, world->r2, tree, &tree[0], 0, world->num);

      /* create threads for parallel tree smoothing length iterators */
      create_smoothing_threads(world, 1, 25, MIN_SMOOTH_LEN, MAX_SMOOTH_LEN, world->r2, tree, &tree[0]);
      //create_smoothing_threads_chunked(world, 1, 25, MIN_SMOOTH_LEN, MAX_SMOOTH_LEN, world->r2, tree, &tree[0]);

      /* create threads for density computation */
      create_density_threads(world);

      /* create threads for pressure computation */
      create_pressure_threads(world);

      /* create threads for sound speed computation */
      create_soundspeed_threads(world);

      /* create threads for CFL computation */
      create_CFL_threads(world);

      /* create threads  for particle time bin update */
      create_timebin_threads(world);

      /* determine minimum timestep to take from maximum bin */
      max_bin=0;
      for(nn=0;nn<world->num;nn++){
	if(world->time_bin[nn]>max_bin)
	  max_bin=world->time_bin[nn];
      }
      world->sub_dt=world->dt/pow(2,max_bin+1);

      /* create threads for hydrodynamic acceleration and internal energy computation */
      create_acceleration_threads(world);
      
      /* timer stop */
      clock_gettime(CLOCK_MONOTONIC, &time2);

      /* compute time */
      timediff(time1, time2, &sph_time);

      /* timer start */
      clock_gettime(CLOCK_MONOTONIC, &time1);

      /* integrate corrector */
      create_corrector_threads(world);

      /* timer stop */
      clock_gettime(CLOCK_MONOTONIC, &time2);

      /* compute time */
      timediff(time1, time2, &int2_time);

      /* compute total, average and minimum energy */
      world->avg_u=0.0;
      world->u_int=0.0;
      world->min_u=HUGE_VAL;
      for(nn=0;nn<world->num;nn++){
	world->u_int+=world->u[nn];
	if(world->u[nn]<world->min_u)
	  world->min_u=world->u[nn];
      }
      world->avg_u=world->u_int/world->num;

      /* compute maximum and minimum of neighbouring particles */
      min_N=world->num;
      max_N=0;
      avg_N=0;
      for(nn=0;nn<world->num;nn++){
	if(world->num_neighbours[nn]>max_N)
	  max_N=world->num_neighbours[nn];

	if(world->num_neighbours[nn]<min_N)
	  min_N=world->num_neighbours[nn];

	avg_N+=world->num_neighbours[nn];
      }
      avg_N=avg_N/world->num;

#if (defined ENERGY_PROFILING)&&ENERGY_PROFILING
      /* compute total gravitational potential and kinetic energy */
      create_total_energy_threads(world, 0.1);

      total_u=total_u2;
      total_u2=world->u_int+world->u_grav+world->u_kin;

      /* display the state of the system */
      printf("time: %.1fyr dt: %f cells: %d avg_N: %.1f u_total: %.1f u_error: %f%%\n",
      	     world->time, world->sub_dt, tree[0].numcells, avg_N, total_u2, 100.0*fabs(total_u2-total_u)/fabs(total_u2));
#else
      /* display current state of the system */
      printf("time: %.1fyr dt: %f cells: %d avg_N: %.1f cputime/yr: %fs\n",
      	     world->time, world->sub_dt, tree[0].numcells, avg_N,
	     ((double)(treetime.tv_sec)+(double)(treetime.tv_nsec)*1.0E-9+
	     (double)(sph_time.tv_sec)+(double)(sph_time.tv_nsec)*1.0E-9+
	     (double)(int_time.tv_sec+int2_time.tv_sec)+(double)(int_time.tv_nsec+int2_time.tv_nsec)*1.0E-9)
	     /(world->sub_dt));
      /*
      printf("time: %.1fyr dt: %f cells: %d avg_N: %.1f tree_t: %fms sph_t: %fms int_t: %fms\n",
      	     world->time, world->sub_dt, tree[0].numcells, avg_N,
	     (double)(treetime.tv_sec)*1.0E3+(double)(treetime.tv_nsec)*1.0E-6,
	     (double)(sph_time.tv_sec)*1.0E3+(double)(sph_time.tv_nsec)*1.0E-6,
	     (double)(int_time.tv_sec+int2_time.tv_sec)*1.0E3+(double)(int_time.tv_nsec+int2_time.tv_nsec)*1.0E-6);
      */
#endif
      /* next time step */
      tt++;

      world->time+=world->sub_dt;
    }

#if ENABLE_GUI
    /* update graphics */
    glClearColor(0.0f,0.0f,0.0f,0.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, 0.0f);
    glRotatef(rot[0]+rot2[0],0.0f,1.0f,0.0f);
    glRotatef(rot[1]+rot2[1],1.0f,0.0f,0.0f);
    glPointSize(1.0f);
    glBegin(GL_POINTS);

    /* find maximum density and energy for normalization */
    max_density=0.0;
    max_pressure=0.0;
    for(nn=0;nn<n;nn++){
      density=world->rho[nn];
      pressure=world->p[nn];

      if(density>max_density)
	max_density=density;

      if(pressure>max_pressure)
	max_pressure=pressure;
    }

    for(nn=0;nn<n;nn++){
      //density=(world->rho[nn]/max_density)*0.8+0.2;
      pressure=(world->p[nn]/max_pressure);
      pressure=pow(pressure,1.0/6.0);
      colormap(pressure,col);
      glColor4f((float)(col->r),(float)(col->g),(float)(col->b),(float)(1.0));
      glVertex3f((zoom)*(world->r[nn*m+0]),
		 (zoom)*(world->r[nn*m+1]),
		 (zoom)*(world->r[nn*m+2]));
    }
    glEnd();
    glFlush();
    SDL_GL_SwapBuffers();

    // write the opengl view as tiff on disk
    //if(fmod(world->time,0.1)<fmod(world->time-world->sub_dt,0.1)){
    //  sprintf(filename, "/usr/crap/testrun/%08d.tif",tiff_frame);
    //  writeframe(filename);
    //  tiff_frame++;
    //}
#endif
  }
    
  printf("\ncleaning up...\n");    
    
  /* clean up vectors */
  free(world->num_neighbours);
  free(a_tree);
  free(a_sph);
  free(r);
  free(world->cellindex);
  free(world->dt_CFL);
  free(world->h);
  free(world->del_rho);
  free(world->u);
  free(world->u2);
  free(world->du);
  free(world->du2);
  free(world->c);
  free(world->p);
  free(world->rho);
  free(world->a2);
  free(world->a);
  free(world->v2);
  free(world->v);
  free(world->r2);
  free(world->r);
  free(world->m);
  free(world->neighbour_list);
  free(world->kick_list);
  free(world->kick);
  free(world->time_bin);

  /* free universe structure */
  free(world);
    
  /* clean up tree */
  free(tree);

#if (defined ENABLE_GUI)&&ENABLE_GUI
  /* clean up SDL */
  SDL_Quit();
#endif

  /* clean up posix threads */
  pthread_exit(NULL);

  return 0;
}

