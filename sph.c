
/* smoothed particle hydrodynamics routines */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dubstep.h"
#include "sph.h"
#include "tree.h"

#ifndef PI
#define PI 3.14159265358979323846264
#endif

double kernel(double q, double h){
  double val=q;

  if(q<1.0){
    val=1.0-3.0/2.0*val*val+3.0/4.0*val*val*val;
  }
  else if(q<2.0){
    val=2.0-val;
    val=1.0/4.0*val*val*val;
  }
  else
    val=0;

  return val/(PI*h*h*h);
}

double kernel_d(double q, double h){
  double val=q;

  if(q<1.0){
    val=-3.0*val+9.0/4.0*val*val;
  }
  else if(q<2.0){
    val=2.0-val;
    val=-3.0/4.0*val*val;
  }
  else
    val=0;

  return val/(PI*h*h*h);
}

double artificial_viscosity(double *dv_ij, double h_ij, double rho_ij,
			    double c_ij, double *dr, double rr, double alpha, double beta, double neta){
  double mu_ij;

  mu_ij=0;
  if((dv_ij[0]*dr[0]+dv_ij[1]*dr[1]+dv_ij[2]*dr[2])<0.0)
    mu_ij=(dv_ij[0]*dr[0]+dv_ij[1]*dr[1]+dv_ij[2]*dr[2])/(rr*rr/(h_ij)+neta);

  return (-alpha*mu_ij*c_ij+beta*mu_ij*mu_ij)/rho_ij;
}

void smooth_velocity_field(struct universe *world, int lo, int hi){
  /* vector norm variables */
  double r;
  double dr[3];

  /* loop variables */
  int ii,jj,kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *rho_in;
  double *m_in;
  double *h_in;
  double *v_in;

  /* smoothed output vector */
  double *v_out;
  
  m=world->dim;  
  n=world->num;

  r_in=world->r2;
  rho_in=world->rho;
  m_in=world->m;
  h_in=world->h;
  v_in=world->v;

  v_out=(double*)malloc(m*n*sizeof(double));
  if(!v_out){
    printf("Out of memory: velocity field smoothing vector not allocated.\n");
    exit(1);
  }

  /* compute filtered velocity field */
  for(ii=lo;ii<hi;ii++){
    v_out[3*ii+0]=0;
    v_out[3*ii+1]=0;
    v_out[3*ii+2]=0;

    n=world->neighbour_list[ii].num;

    for(jj=0;jj<n;jj++){
      kk=world->neighbour_list[ii].list[jj];

      /* particle-particle distance */
      dr[0]=r_in[3*ii+0]-r_in[3*kk+0];
      dr[1]=r_in[3*ii+1]-r_in[3*kk+1];
      dr[2]=r_in[3*ii+2]-r_in[3*kk+2];
      r=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])/h_in[ii];

      v_out[3*ii+0]+=1/rho_in[ii]*m_in[kk]*v_in[3*kk+0]*0.5*(kernel(r,h_in[ii])+kernel(r,h_in[kk]));
      v_out[3*ii+1]+=1/rho_in[ii]*m_in[kk]*v_in[3*kk+1]*0.5*(kernel(r,h_in[ii])+kernel(r,h_in[kk]));
      v_out[3*ii+2]+=1/rho_in[ii]*m_in[kk]*v_in[3*kk+2]*0.5*(kernel(r,h_in[ii])+kernel(r,h_in[kk]));
    }
  }

  /* free original velocity vector and replace with filtered */
  free(world->v);
  world->v=v_out;
}

void smooth_energy_field(struct universe *world, int lo, int hi){
  /* vector norm variables */
  double r;
  double dr[3];

  /* loop variables */
  int ii,jj,kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *rho_in;
  double *m_in;
  double *h_in;
  double *u_in;

  /* smoothed output vector */
  double *u_out;
  
  m=world->dim;  
  n=world->num;

  r_in=world->r2;
  rho_in=world->rho;
  m_in=world->m;
  h_in=world->h;
  u_in=world->u;

  u_out=(double*)malloc(n*sizeof(double));
  if(!u_out){
    printf("Out of memory: energy smoothing vector not allocated.\n");
    exit(1);
  }

  /* compute filtered energy field */
  for(ii=lo;ii<hi;ii++){
    u_out[ii]=0;

    n=world->neighbour_list[ii].num;

    for(jj=0;jj<n;jj++){
      kk=world->neighbour_list[ii].list[jj];

      /* particle-particle distance */
      dr[0]=r_in[3*ii+0]-r_in[3*kk+0];
      dr[1]=r_in[3*ii+1]-r_in[3*kk+1];
      dr[2]=r_in[3*ii+2]-r_in[3*kk+2];
      r=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])/h_in[ii];

      u_out[ii]+=1/rho_in[ii]*m_in[kk]*u_in[kk]*0.5*(kernel(r,h_in[ii])+kernel(r,h_in[kk]));
    }
  }

  /* free original velocity vector and replace with filtered */
  free(world->u);
  world->u=u_out;
}

void compute_smoothing_length_tree(struct universe *world, double min_h, double max_h, int iterations, int N_target,
				   double *r, struct cell *tree, struct cell *root, int lo, int hi){
  /* loop variables */
  int ii,jj;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *h_in;

  /* pointer to particle neighbour number vector */
  int *num_neighbours_in;

  /* pointer to smoothing length iteration variables */
  double h,h_new;

  /* number of interacting particles */
  int N;

  /* particle list buffer */
  int *buffer;

  m=world->dim;  
  n=world->num;

  r_in=r;
  h_in=world->h;

  num_neighbours_in=world->num_neighbours;

  buffer=(int*)malloc(n*sizeof(int));
  if(!buffer){
    printf("Out of memory: smoothing length iteration buffer not allocated.\n");
    exit(1);
  }

  /* iterate towards optimum number of neighbours */
  for(ii=lo;ii<hi;ii++){
    h=h_in[ii];
      
    for(jj=0;jj<iterations;jj++){
      N=0;
      neighbour_walk(tree, root, &r_in[3*ii], h, max_h, h_in, &N, buffer);
      
      h_new=h*0.5*(1.0+pow((double)(N_target)/(double)(N),1.0/3.0));
      if(h_new>min_h&&h_new<max_h)
	h=h_new;
      else
	break;
    }
    h_in[ii]=h;
    num_neighbours_in[ii]=N;
    
    world->neighbour_list[ii].num=N;
    if(world->neighbour_list[ii].list)
      free(world->neighbour_list[ii].list);
    
    world->neighbour_list[ii].list=(int*)malloc(N*sizeof(int));
    if(!world->neighbour_list[ii].list){
      printf("Out of memory: particle neighbour list not allocated.\n");
      exit(1);
    }
    memcpy(world->neighbour_list[ii].list, buffer, N*sizeof(int));
  }
  free(buffer);
}

void compute_constant_smoothing_length_tree(struct universe *world, double min_h, double max_h, int iterations, int N_target,
					    double *r, struct cell *tree, struct cell *root, int lo, int hi){
  /* loop variables */
  int ii;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *h_in;

  /* pointer to particle neighbour number vector */
  int *num_neighbours_in;

  /* number of interacting particles */
  int N;

  /* particle list buffer */
  int *buffer;

  m=world->dim;  
  n=world->num;

  r_in=r;
  h_in=world->h;

  num_neighbours_in=world->num_neighbours;

  buffer=(int*)malloc(n*sizeof(int));
  if(!buffer){
    printf("Out of memory: smoothing length iteration buffer not allocated.\n");
    exit(1);
  }

  /* iterate towards optimum number of neighbours */
  for(ii=lo;ii<hi;ii++){
    N=0;
    neighbour_walk(tree, root, &r_in[3*ii], h_in[ii], max_h, h_in, &N, buffer);
      
    num_neighbours_in[ii]=N;
    
    world->neighbour_list[ii].num=N;
    if(world->neighbour_list[ii].list)
      free(world->neighbour_list[ii].list);
    
    world->neighbour_list[ii].list=(int*)malloc(N*sizeof(int));
    if(!world->neighbour_list[ii].list){
      printf("Out of memory: particle neighbour list not allocated.\n");
      exit(1);
    }
    memcpy(world->neighbour_list[ii].list, buffer, N*sizeof(int));
  }
  free(buffer);
}

void compute_smoothing_length_neighbours(struct universe *world, int iterations, int N_target, int lo, int hi){
  /* vector norm variables */
  double r;
  double dr[3];

  /* loop variables */
  int ii,jj,kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *h_in;

  int *num_neighbours_in;

  double h,h_new;

  /* number of interacting particles */
  int N;

  /* particle list buffer */
  int *buffer;

  m=world->dim;  
  n=world->num;

  r_in=world->r2;
  h_in=world->h;

  num_neighbours_in=world->num_neighbours;

  buffer=(int*)malloc(2*n*sizeof(int));
  if(!buffer){
    printf("Out of memory: smoothing length iteration buffer not allocated.\n");
    exit(1);
  }

  /* iterate towards optimum number of neighbours */
  for(ii=lo;ii<hi;ii++){
    h=h_in[ii];
      
    if(world->neighbour_list[ii].list){
      world->neighbour_list[ii].num=0;
      free(world->neighbour_list[ii].list);
    }

    for(kk=0;kk<iterations;kk++){
      N=0;
      for(jj=0;jj<n;jj++){
	/* particle-particle distance */
	dr[0]=r_in[3*ii+0]-r_in[3*jj+0];
	dr[1]=r_in[3*ii+1]-r_in[3*jj+1];
	dr[2]=r_in[3*ii+2]-r_in[3*jj+2];
	r=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	if(r/h<2.0||r/h_in[jj]<2.0){
	  buffer[N]=jj;
	  N++;
	}
      }
      h_new=h*0.5*(1+pow((double)(N_target)/(double)(N),1.0/3.0));
      if(h_new>0.01&&h_new<1.0)
	h=h_new;
      else
	break;
    }
    h_in[ii]=h;
    num_neighbours_in[ii]=N;
    
    world->neighbour_list[ii].num=N;
    world->neighbour_list[ii].list=(int*)malloc(N*sizeof(int));
    if(!world->neighbour_list[ii].list){
      printf("Out of memory: particle neighbour list not allocated.\n");
      exit(1);
    }
    memcpy(world->neighbour_list[ii].list, buffer, N*sizeof(int));
  }
  free(buffer);
}

void compute_smoothing_length_mass(struct universe *world, int iterations, double m_target, int lo, int hi){
  /* vector norm variables */
  double r;
  double dr[3];

  /* loop variables */
  int ii,jj,kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *h_in;
  double *m_in;

  double *dt_CFL_in;

  int *num_neighbours_in;

  double h,h_new;

  /* neighboring mass */
  double N_mass;

  int nn;

  /* particle list buffer */
  int *buffer;

  m=world->dim;  
  n=world->num;

  r_in=world->r2;
  h_in=world->h;
  m_in=world->m;
  dt_CFL_in=world->dt_CFL;

  num_neighbours_in=world->num_neighbours;

  buffer=(int*)malloc(2*n*sizeof(int));

  /* iterate towards optimum number of neighbours */
  for(ii=lo;ii<hi;ii++){
    h=h_in[ii];

    if(world->neighbour_list[ii].list){
      world->neighbour_list[ii].num=0;
      free(world->neighbour_list[ii].list);
    }

    for(kk=0;kk<iterations;kk++){
      N_mass=0;
      nn=0;
      for(jj=0;jj<n;jj++){
	/* particle-particle distance */
	dr[0]=r_in[3*ii+0]-r_in[3*jj+0];
	dr[1]=r_in[3*ii+1]-r_in[3*jj+1];
	dr[2]=r_in[3*ii+2]-r_in[3*jj+2];
	r=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

	if(r/h<2.0||r/h_in[jj]<2.0){
	  buffer[nn++]=jj;
	  N_mass+=m_in[jj];
	}
      }
      h_new=h*0.5*(1+pow((double)(m_target)/(double)(N_mass),1.0/3.0));
      if(h_new>0.01&&h_new<4.0)
	h=h_new;
      else
	break;
    }
    h_in[ii]=h;
    num_neighbours_in[ii]=nn;
    
    world->neighbour_list[ii].num=nn;
    world->neighbour_list[ii].list=(int*)malloc(nn*sizeof(int));
    if(!world->neighbour_list[ii].list){
      printf("Out of memory: particle neighbour list not allocated.\n");
      exit(1);
    }
    memcpy(world->neighbour_list[ii].list, buffer, nn*sizeof(int));
  }
  free(buffer);
}

void compute_density(struct universe *world, double h, int lo, int hi){
  /* vector norm variables */
  double r;
  double dr[3];

  /* loop variables */
  int ii,jj,kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *rho_in;
  double *m_in;
  double *h_in;
  
  m=world->dim;  
  n=world->num;

  r_in=world->r2;
  rho_in=world->rho;
  m_in=world->m;
  h_in=world->h;

  /* compute density */
  for(ii=lo;ii<hi;ii++){
    rho_in[ii]=0;
    n=world->neighbour_list[ii].num;
    for(jj=0;jj<n;jj++){
      kk=world->neighbour_list[ii].list[jj];
      /* particle-particle distance */
      dr[0]=r_in[3*ii+0]-r_in[3*kk+0];
      dr[1]=r_in[3*ii+1]-r_in[3*kk+1];
      dr[2]=r_in[3*ii+2]-r_in[3*kk+2];
      r=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2])/h_in[ii];
      rho_in[ii]+=m_in[kk]*kernel(r,h_in[ii]);
    }
  }
}

void compute_pressure(struct universe *world, double gamma, int lo, int hi){
  /* loop variables */
  int ii;

  /* state vector dimensions */
  int n;

  /* pointers to state vectors */
  double *p_in;
  double *rho_in;
  double *u_in;

  n=world->num;

  p_in=world->p;
  rho_in=world->rho;
  u_in=world->u2;

  /* compute pressure using polytropic equation of state */
  for(ii=lo;ii<hi;ii++){
    p_in[ii]=(gamma-1)*rho_in[ii]*u_in[ii];
  }
}

void compute_soundspeed(struct universe *world, double gamma, int lo, int hi){
  /* loop variables */
  int ii;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *p_in;
  double *rho_in;
  double *c_in;

  m=world->dim;  
  n=world->num;

  p_in=world->p;
  rho_in=world->rho;
  c_in=world->c;

  for(ii=lo;ii<hi;ii++){
    c_in[ii]=sqrt(gamma*p_in[ii]/rho_in[ii]);
  }
}

void compute_cfl(struct universe *world, double C_0, int lo, int hi){
  /* vector norm variables */
  double rr;
  double dr[3];

  /* loop variables */
  int ii,jj,kk;

  /* smoothing kernel variable */
  double dW;
  double dW_ij[3];

  double dv_ij[3];

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *rho_in;
  double *m_in;
  double *h_in;
  double *v_in;
  double *c_in;

  double *dt_CFL_in;

  double div_v_ij;

  double mu_ij,max_mu;

  double h_ij;

  double alpha;
  double beta;
  double neta=0.01;
  
  m=world->dim;  
  n=world->num;

  r_in=world->r2;
  rho_in=world->rho;
  m_in=world->m;
  h_in=world->h;
  v_in=world->v2;
  c_in=world->c;
  dt_CFL_in=world->dt_CFL;

  alpha=world->alpha;
  beta=world->beta;

  for(ii=lo;ii<hi;ii++){
    div_v_ij=0;
    max_mu=0;
    n=world->neighbour_list[ii].num;
    
    for(jj=0;jj<n;jj++){
      kk=world->neighbour_list[ii].list[jj];
      
      if(ii!=kk){
	/* particle-particle distance */
	dr[0]=r_in[3*ii+0]-r_in[3*kk+0];
	dr[1]=r_in[3*ii+1]-r_in[3*kk+1];
	dr[2]=r_in[3*ii+2]-r_in[3*kk+2];
	rr=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	
	dv_ij[0]=v_in[3*ii+0]-v_in[3*kk+0];
	dv_ij[1]=v_in[3*ii+1]-v_in[3*kk+1];
	dv_ij[2]=v_in[3*ii+2]-v_in[3*kk+2];
	
	dW=0.5*(kernel_d(rr/h_in[ii],h_in[ii])+kernel_d(rr/h_in[kk],h_in[kk]));
	
	dW_ij[0]=dW*dr[0]/rr;
	dW_ij[1]=dW*dr[1]/rr;
	dW_ij[2]=dW*dr[2]/rr;
	
	div_v_ij-=m_in[kk]/rho_in[ii]*(dv_ij[0]*dW_ij[0]+dv_ij[1]*dW_ij[1]+dv_ij[2]*dW_ij[2]);
	
	h_ij=0.5*(h_in[ii]+h_in[kk]);
	
	mu_ij=0;
	if((dv_ij[0]*dr[0]+dv_ij[1]*dr[1]+dv_ij[2]*dr[2])<0)
	  mu_ij=(dv_ij[0]*dr[0]+dv_ij[1]*dr[1]+dv_ij[2]*dr[2])/(rr*rr/(h_ij)+neta);
	
	if(abs(mu_ij)>max_mu)
	  max_mu=abs(mu_ij);
      }
    }
    dt_CFL_in[ii]=(C_0*h_in[ii])/(h_in[ii]*abs(div_v_ij)+c_in[ii]+1.2*(alpha*c_in[ii]+beta*max_mu));
  }
}

void compute_internal_energy(struct universe *world, double *a_sph, int lo, int hi){
  /* vector norm variables */
  double rr;
  double dr[3];

  /* smoothing kernel variable */
  double dW;
  double dW_ij[3];

  double dv_ij[3];

  double t;

  /* artificial viscosity variables */
  double Pi_ij;
  double alpha=1.0;
  double beta=1.0;
  double neta=0.01;

  /* loop variables */
  int ii;
  int jj;
  int kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *v_in;
  double *rho_in;
  double *c_in;
  double *du_in;
  double *m_in;
  double *p_in;
  double *h_in;
  double *dt_CFL_in;

  m=world->dim;  
  n=world->num;

  r_in=world->r2;
  v_in=world->v2;
  rho_in=world->rho;
  c_in=world->c;
  du_in=world->du;
  m_in=world->m;
  p_in=world->p;
  h_in=world->h;
  dt_CFL_in=world->dt_CFL;

  alpha=world->alpha;
  beta=world->beta;

  for(ii=lo;ii<hi;ii++){
    /* init dot energy*/
    du_in[ii]=0;

    n=world->neighbour_list[ii].num;

    for(jj=0;jj<n;jj++){
      kk=world->neighbour_list[ii].list[jj];

      /* compute gradient of kernel */
      if(ii!=kk){
	/* particle-particle distance */
	dr[0]=r_in[3*ii+0]-r_in[3*kk+0];
	dr[1]=r_in[3*ii+1]-r_in[3*kk+1];
	dr[2]=r_in[3*ii+2]-r_in[3*kk+2];
	rr=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
      
	dv_ij[0]=v_in[3*ii+0]-v_in[3*kk+0];
	dv_ij[1]=v_in[3*ii+1]-v_in[3*kk+1];
	dv_ij[2]=v_in[3*ii+2]-v_in[3*kk+2];
	
	/* compute artificial viscosity term */
	Pi_ij=artificial_viscosity(dv_ij,
				   0.5*(h_in[ii]+h_in[kk]),
				   0.5*(rho_in[ii]+rho_in[kk]),
				   0.5*(c_in[ii]+c_in[kk]),
				   dr, rr,
				   alpha, beta, neta);
	
	/* asymmetric */
	t=m_in[kk]*(p_in[ii]/(rho_in[ii]*rho_in[ii])+0.5*Pi_ij);
	
	dW=0.5*(kernel_d(rr/h_in[ii],h_in[ii])+kernel_d(rr/h_in[kk],h_in[kk]));
	
	dW_ij[0]=dW*dr[0]/rr;
	dW_ij[1]=dW*dr[1]/rr;
	dW_ij[2]=dW*dr[2]/rr;
	
	du_in[ii]+=t*(dv_ij[0]*dW_ij[0]+dv_ij[1]*dW_ij[1]+dv_ij[2]*dW_ij[2]);
      }
    }
  }
}

void compute_sph_acceleration(struct universe *world, double *r, double *v, double *a, int lo, int hi){
  /* vector norm variables */
  double rr;
  double dr[3];

  /* smoothing kernel variable */
  double dW;
  double dW_ij[3];

  double dv_ij[3];

  double t;

  /* artificial viscosity variables */
  double Pi_ij;
  double alpha=1.0;
  double beta=1.0;
  double neta=0.01;

  /* loop variables */
  int ii;
  int jj;
  int kk;

  /* state vector dimensions */
  int m;
  int n;

  /* pointers to state vectors */
  double *r_in;
  double *v_in;
  double *rho_in;
  double *c_in;
  double *du_in;
  double *m_in;
  double *p_in;
  double *h_in;
  double *dt_CFL_in;

  m=world->dim;  
  n=world->num;

  r_in=r;
  v_in=v;

  rho_in=world->rho;
  c_in=world->c;
  du_in=world->du;
  m_in=world->m;
  p_in=world->p;
  h_in=world->h;
  dt_CFL_in=world->dt_CFL;

  alpha=world->alpha;
  beta=world->beta;

  for(ii=lo;ii<hi;ii++){
    if(world->kick[ii]){
      /* init acceleration vector */
      a[3*ii+0]=0;
      a[3*ii+1]=0;
      a[3*ii+2]=0;
    
      n=world->neighbour_list[ii].num;
    
      for(jj=0;jj<n;jj++){
	kk=world->neighbour_list[ii].list[jj];
      
	/* compute gradient of kernel */
	if(ii!=kk){
	  /* particle-particle distance */
	  dr[0]=r_in[3*ii+0]-r_in[3*kk+0];
	  dr[1]=r_in[3*ii+1]-r_in[3*kk+1];
	  dr[2]=r_in[3*ii+2]-r_in[3*kk+2];
	  rr=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
	  
	  dv_ij[0]=v_in[3*ii+0]-v_in[3*kk+0];
	  dv_ij[1]=v_in[3*ii+1]-v_in[3*kk+1];
	  dv_ij[2]=v_in[3*ii+2]-v_in[3*kk+2];
	
	  /* compute artificial viscosity term */
	  Pi_ij=artificial_viscosity(dv_ij,
				     0.5*(h_in[ii]+h_in[kk]),
				     0.5*(rho_in[ii]+rho_in[kk]),
				     0.5*(c_in[ii]+c_in[kk]),
				     dr, rr,
				     alpha, beta, neta);
	
	  /* geometric mean symmetrization */
	  //t=m_in[kk]*(2.0*sqrt(p_in[ii]*p_in[kk])/(rho_in[ii]*rho_in[kk])+Pi_ij);
	
	  /* arithmetic mean symmetrization */
	  t=m_in[kk]*(p_in[ii]/(rho_in[ii]*rho_in[ii])+p_in[kk]/(rho_in[kk]*rho_in[kk])+Pi_ij);
	
	  dW=0.5*(kernel_d(rr/h_in[ii],h_in[ii])+kernel_d(rr/h_in[kk],h_in[kk]));
	  
	  dW_ij[0]=dW*dr[0]/rr;
	  dW_ij[1]=dW*dr[1]/rr;
	  dW_ij[2]=dW*dr[2]/rr;
	  
	  a[3*ii+0]-=t*dW_ij[0];
	  a[3*ii+1]-=t*dW_ij[1];
	  a[3*ii+2]-=t*dW_ij[2];
	}
      }
    }
  }
}

