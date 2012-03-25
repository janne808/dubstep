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
    TIFFSetField(file, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
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

  /* initial timestep, determined from CFL before main loop */
  double dt=0.0;

  /* plummer gravitational softening factor*/
  double epsilon=0.005;

  /* initial smoothing length */
  double h=1.0;

  /* artificial viscosity parameters */
  double alpha=1.0;
  double beta=1.0;

  /* adiabatic exponent*/
  double gamma=1.4;

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

  double energy;
  double max_energy;

  struct color *col=0;
#endif

  double avg_u;
  double min_u;
  double min_dt;

#if ENABLE_GUI
  char filename[128];
  double wtime=0;
  int wtt=0;
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
  world->last_kick=(double*)malloc(n*sizeof(double));

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
    world->last_kick[ii]=0.0;
    world->m[ii]=1.0/(float)(n);
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
    world->u[ii]=4.0;
    world->u2[ii]=world->u[ii];
  }

  /* initial displacement */
  for(ii=0;ii<n;ii++){
    double x,y,z,rr;

    x=2*((double)(rand())/RAND_MAX)-1;
    y=2*((double)(rand())/RAND_MAX)-1;
    z=2*((double)(rand())/RAND_MAX)-1;

    while(sqrt(x*x+y*y+z*z)>1.0){
      x=2*((double)(rand())/RAND_MAX)-1;
      y=2*((double)(rand())/RAND_MAX)-1;
      z=2*((double)(rand())/RAND_MAX)-1;
    }

    world->r[ii*m+0]=24.0*x;
    world->r[ii*m+1]=24.0*y;
    world->r[ii*m+2]=24.0*z;
    world->r2[ii*m+0]=world->r[ii*m+0];
    world->r2[ii*m+1]=world->r[ii*m+1];
    world->r2[ii*m+2]=world->r[ii*m+2];

    rr=sqrt(x*x+y*y+z*z);

    world->v[ii*m+0]=-1.0*y/rr;
    world->v[ii*m+1]=1.0*x/rr;
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
  createSmoothingThreads(world, 10, 25);

  /* create threads for density computation */
  createDensityThreads(world);

  /* create threads for pressure computation */
  createPressureThreads(world);

  /* create threads for sound speed computation */
  createSoundspeedThreads(world);

  /* create threads for CFL computation */
  createCFLThreads(world);

  /* determine minimum timestep to take from CFL */
  min_dt=HUGE_VAL;
  for(nn=0;nn<world->num;nn++){
    if(world->dt_CFL[nn]<min_dt)
      min_dt=world->dt_CFL[nn];
  }
  world->sub_dt=min_dt;
  
  /* create threads for sph energy computation */
  createEnergyThreads(world);

  /* create threads for hydrodynamic acceleration computation */
  createAccelerationThreads(world);
      
  /* initialize tree root parameters */
  init_treeroot(tree, world, world->r);
    
  /* form tree */
  branchrecurse(tree, &tree[0], world->cellindex);
      
  /* integrate corrector */
  createCorrectorThreads(world);

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
      createPredictorThreads(world);

      t2=SDL_GetTicks();
      int_time=t2-t1;

      /* compute sph variables */

      /* create threads for smoothing length interation and
	 interacting particle list generation */
      //createSmoothingThreads(world, 1, 25);

      t1=SDL_GetTicks();

      /* initialize tree root parameters */
      init_treeroot(tree, world, world->r2);
    
      /* form tree */
      branchrecurse(tree, &tree[0], world->cellindex);

      /* serial tree smoothing length iterator */
      compute_smoothing_length_tree(world, h, 1, 25, world->r2, tree, &tree[0]);

      t2=SDL_GetTicks();
      treetime=t2-t1;

      /* create threads for density computation */
      createDensityThreads(world);

      /* create threads for pressure computation */
      createPressureThreads(world);

      /* create threads for sound speed computation */
      createSoundspeedThreads(world);

      /* create threads for CFL computation */
      createCFLThreads(world);

      /* determine minimum timestep to take from CFL */
      min_dt=HUGE_VAL;
      for(nn=0;nn<world->num;nn++){
	if(world->dt_CFL[nn]<min_dt)
	  min_dt=world->dt_CFL[nn];
      }
      world->sub_dt=min_dt;

      /* create threads for sph energy computation */
      createEnergyThreads(world);

      /* create threads for hydrodynamic acceleration computation */
      createAccelerationThreads(world);
      
      t2=SDL_GetTicks();
      sph_time=t2-t1;

      t1=SDL_GetTicks();

      /* integrate corrector */
      createCorrectorThreads(world);

      t2=SDL_GetTicks();
      int_time+=t2-t1;

      /* compute average and minimum energy */
      avg_u=0.0;
      min_u=HUGE_VAL;
      for(nn=0;nn<world->num;nn++){
	avg_u+=world->u[nn];
	if(world->u[nn]<min_u)
	  min_u=world->u[nn];
      }
      avg_u=avg_u/world->num;

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

      printf("time: %f dt: %f cells: %d tree t: %d ms sph t: %d ms integration t: %d ms avg_u: %f min_u: %f avg_N: %f\n",
      	     world->time, world->sub_dt, tree[0].numcells, treetime, sph_time, int_time, avg_u, min_u, avg_N);

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
    max_energy=0.0;
    for(nn=0;nn<n;nn++){
      density=world->rho[nn];
      energy=world->u[nn];

      if(density>max_density)
	max_density=density;

      if(energy>max_energy)
	max_energy=energy;
    }

    for(nn=0;nn<n;nn++){
      density=(world->rho[nn]/max_density)*0.8+0.2;
      energy=(world->u[nn]/max_energy);
      colormap(energy,col);
      glColor4f((float)(col->r),(float)(col->g),(float)(col->b),(float)(density));
      glVertex3f((zoom)*(world->r[nn*m+0]),
		 (zoom)*(world->r[nn*m+1]),
		 (zoom)*(world->r[nn*m+2]));
    }
    glEnd();
    glFlush();
    SDL_GL_SwapBuffers();

    // write the opengl view as tiff on disk
    //wtime+=world->sub_dt;
    //if(wtime>0.01){
    //  sprintf(filename, "/home/janne808/testrun/%08d.tif",wtt);
    //  writeframe(filename);
    //  while(wtime>0.01)
    //	wtime-=0.01;
    //  wtt++;
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
  free(world->last_kick);
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

