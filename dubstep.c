/* dubstep routines */

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
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include <SDL.h>

#if ENABLE_GUI
#include <SDL_opengl.h>
#include <tiffio.h>
#endif

#include "dubstep.h"
#include "tree.h"
#include "sph.h"
#include "threads.h"

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
  double theta=0.5;

  /* maximum timestep */
  double dt=0.1;

  /* plummer gravitational softening factor*/
  double epsilon=1.0;

  /* initial smoothing length */
  double h=1.0;

  /* artificial viscosity parameters */
  double alpha=1.0;
  double beta=1.0;

  /* adiabatic exponent*/
  double gamma=1.4;

  /* time bin handling variables */
  int max_bin;
  double ksi;
  double old_dt;

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
  unsigned int t1,t2;
  unsigned int treetime;
  unsigned int int_time;
  unsigned int sph_time;

  /* particle neighbour statistics */
  int max_N;
  int min_N;
  double avg_N;

#if ENABLE_GUI
  /* UI variables */
  double rot[3];
  double rot2[3];
  int mouse_x,mouse_y;
  int mouse_dx,mouse_dy;
  int buttonstate;
  double zoom=1/12.0;
  double zoom2;

  double density;
  double max_density;

  double pressure;
  double max_pressure;

  struct color *col=0;
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
    world->neighbour_list[ii].list=0;
    world->neighbour_list[ii].num=0;
    world->num_neighbours[ii]=0;
    world->cellindex[ii]=0;
    world->h[ii]=h;
    world->dt_CFL[ii]=world->sub_dt;
    world->kick[ii]=1;
    world->time_bin[ii]=0;
    world->m[ii]=5.0/(float)(n);
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
    world->u[ii]=1.0;
    world->u2[ii]=world->u[ii];
  }

  /* initial displacement */
  for(ii=0;ii<n;ii++){
    double x,y,z,rr;
    double vv;

    /* box-muller transform for normally distributed samples */
    x=sqrt(-2*log((double)rand()/RAND_MAX))*cos(2.0*PI*(double)rand()/RAND_MAX);
    y=sqrt(-2*log((double)rand()/RAND_MAX))*cos(2.0*PI*(double)rand()/RAND_MAX);
    z=sqrt(-2*log((double)rand()/RAND_MAX))*cos(2.0*PI*(double)rand()/RAND_MAX);

    world->r[ii*m+0]=10.0*x;
    world->r[ii*m+1]=10.0*y;
    world->r[ii*m+2]=2.0*z;
    world->r2[ii*m+0]=world->r[ii*m+0];
    world->r2[ii*m+1]=world->r[ii*m+1];
    world->r2[ii*m+2]=world->r[ii*m+2];

    rr=sqrt(x*x+y*y+z*z);

    vv=sqrt(2.0*G*0.1*2.0/rr/2.0);

    world->v[ii*m+0]=-vv*y/rr;
    world->v[ii*m+1]=vv*x/rr;
    world->v[ii*m+2]=0;    
    //world->v[ii*m+0]=0;
    //world->v[ii*m+1]=0;
    //world->v[ii*m+2]=0;    
    world->v2[ii*m+0]=world->v[ii*m+0];
    world->v2[ii*m+1]=world->v[ii*m+1];
    world->v2[ii*m+2]=world->v[ii*m+2];    
  }

  /* compute initial sph variables */

  /* create threads for smoothing length interation and
     interacting particle list generation */
  //create_smoothing_threads(world, 10, 25);
      
  /* initialize tree root parameters */
  init_treeroot(tree, world, world->r2);
    
  /* form tree */
  branch_recurse(tree, &tree[0], world->cellindex);
      
  /* serial tree smoothing length iterator */
  compute_smoothing_length_tree(world, 1.0, 2.0, 10, 25, world->r2, tree, &tree[0], 0, world->num);

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
  
  /* create threads for sph energy computation */
  create_energy_threads(world);

  /* filter initial energy field */
  smooth_energy_field(world, 0, world->num);

  /* create threads for hydrodynamic acceleration computation */
  create_acceleration_threads(world);

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
      t1=SDL_GetTicks();

      /* integrate predictor */
      create_predictor_threads(world);

      t2=SDL_GetTicks();
      int_time=t2-t1;

      /* compute sph variables */

      /* create threads for smoothing length interation and
	 interacting particle list generation */

      t1=SDL_GetTicks();

      /* initialize tree root parameters */
      init_treeroot(tree, world, world->r2);
    
      /* form tree */
      branch_recurse(tree, &tree[0], world->cellindex);

      /* serial tree smoothing length iterator */
      //compute_smoothing_length_tree(world, h, 1, 25, world->r2, tree, &tree[0], 0, world->num);

      /* create threads for parallel tree smoothing length iterators */
      create_smoothing_threads(world, 1, 10, 1.0, 2.0, world->r2, tree, &tree[0]);

      t2=SDL_GetTicks();
      treetime=t2-t1;

      /* create threads for density computation */
      create_density_threads(world);

      /* create threads for pressure computation */
      create_pressure_threads(world);

      /* create threads for sound speed computation */
      create_soundspeed_threads(world);

      /* create threads for CFL computation */
      create_CFL_threads(world);

      /* update particle time bins */
      for(nn=0;nn<world->num;nn++){
	/* preserve current time step */
	old_dt=world->dt/pow(2,world->time_bin[nn]);

	for(ii=0;ii<world->num;ii++){
	  if(world->dt_CFL[nn]>world->dt/pow(2,ii)){
	    /* particle can always move to a larger bin */
	    if(world->time_bin[nn]<ii){
	      //printf("Particle %d moves from bin %d to %d.\n", nn, world->time_bin[nn], ii);
	      ksi=(world->dt/pow(2,world->time_bin[nn]))/(world->dt/pow(2,ii));
	      world->time_bin[nn]=ii;
	    }
	    /* move to smaller bin only if the current bin is in sync with it */
	    else if(world->time_bin[nn]>ii){
	      if(floor(world->time/(world->dt/pow(2,world->time_bin[nn])))>
		 floor((world->time-world->sub_dt)/(world->dt/pow(2,world->time_bin[nn])))&&
		 floor(world->time/(world->dt/pow(2,world->time_bin[nn]+1)))>
		 floor((world->time-world->sub_dt)/(world->dt/pow(2,world->time_bin[nn]+1)))
		 ){
		//printf("Particle %d moves from bin %d to %d.\n", nn, world->time_bin[nn], ii);
		ksi=(world->dt/pow(2,world->time_bin[nn]))/(world->dt/pow(2,ii));
		world->time_bin[nn]=ii;
	      }
	    }
	    /* if time bin doesnt change, break before position correction */
	    else
	      break;

	    /* compute position correction */
	    world->r2[3*nn+0]-=(1-1/ksi)*(1+1/ksi)*old_dt*old_dt/8.0*world->a2[3*nn+0];
	    world->r2[3*nn+1]-=(1-1/ksi)*(1+1/ksi)*old_dt*old_dt/8.0*world->a2[3*nn+1];
	    world->r2[3*nn+2]-=(1-1/ksi)*(1+1/ksi)*old_dt*old_dt/8.0*world->a2[3*nn+2];
	    break;
	  }
	}
	//printf("%d\t", world->time_bin[nn]);
      }
      //printf("\n");

      /* determine which particles need to be kicked with the old sub_dt value */
      for(nn=0;nn<world->num;nn++){
	if(floor(world->time/(world->dt/pow(2,world->time_bin[nn])))>
	   floor((world->time-world->sub_dt)/(world->dt/pow(2,world->time_bin[nn]))))
	  world->kick[nn]=1;
	else
	  world->kick[nn]=0;
	//printf("%d\t", world->kick[nn]);	
      }
      //printf("\n");

      /* determine minimum timestep to take from maximum bin */
      max_bin=0;
      for(nn=0;nn<world->num;nn++){
	if(world->time_bin[nn]>max_bin)
	  max_bin=world->time_bin[nn];
      }
      world->sub_dt=world->dt/pow(2,max_bin+1);

      /* create threads for sph energy computation */
      create_energy_threads(world);

      /* create threads for hydrodynamic acceleration computation */
      create_acceleration_threads(world);
      
      t2=SDL_GetTicks();
      sph_time=t2-t1;

      t1=SDL_GetTicks();

      /* integrate corrector */
      create_corrector_threads(world);

      t2=SDL_GetTicks();
      int_time+=t2-t1;

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

      /* compute total gravitational potential and kinetic energy */
      //create_total_energy_threads(world, 0.1);

      /* display the state of the system */
      printf("time: %f dt: %f cells: %d avg_N: %f u_int: %f u_grav: %f u_kin: %f\n",
      	     world->time, world->sub_dt, tree[0].numcells, avg_N, world->u_int, world->u_grav, world->u_kin);

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
      density=(world->rho[nn]/max_density)*0.8+0.2;
      pressure=(world->p[nn]/max_pressure);
      pressure=pow(pressure,1.0/3.0);
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
    //sprintf(filename, "/usr/crap/testrun/%08d.tif",tt);
    //writeframe(filename);
#endif
  }
    
  printf("\ncleaning up...\n");    
    
  /* clean up vectors */
  free(world->num_neighbours);
  free(a_tree);
  free(a_sph);
  free(r);
  free(world->cellindex);
  free(world->kick);
  free(world->time_bin);
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

  /* free universe structure */
  free(world);
    
  /* clean up tree */
  free(tree);

  /* clean up SDL */
  SDL_Quit();

  /* clean up posix threads */
  pthread_exit(NULL);

  return 0;
}

