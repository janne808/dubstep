/* barnes-hut tree routines */

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
#include "sph.h"
#include "tree.h"
#include "linear.h"

int init_treeroot(struct cell *tree, struct universe *world, dubfloat_t *r){
  int ii;
  int n,m;
  
  /* determine state vector dimensions */
  n=world->num;
  m=world->dim;

  /* form tree root cell*/
  /* init particle displacement vector */
  tree[0].r=malloc(m*n*sizeof(dubfloat_t));
  if(!tree[0].r){
    printf("Out of memory: root displacement vector not allocated.\n");
    return -1;
  }
  memcpy(tree[0].r, r, m*n*sizeof(dubfloat_t));

  /* init particle mass vector */
  tree[0].m=malloc(n*sizeof(dubfloat_t));        
  if(!tree[0].m){
    printf("Out of memory: root mass vector not allocated.\n");
    return -1;
  }
  memcpy(tree[0].m, world->m , n*sizeof(dubfloat_t));
  
  /* init particle index vector */
  tree[0].particle_index_list=malloc(n*sizeof(int));
  if(!tree[0].particle_index_list){
    printf("Out of memory: root particle index vector not allocated.\n");
    return -1;
  }
  for(ii=0;ii<n;ii++){
    tree[0].particle_index_list[ii]=ii;
  }

  /* init particle cell index vector */
  for(ii=0;ii<n;ii++){
    world->cellindex[ii]=-1;
  }

  /* search for maximum values */
  /* make sure we are dealing with squares */
  tree[0].space[m*0+0]=-INFINITY;
  tree[0].space[m*0+1]=-INFINITY;
  tree[0].space[m*0+2]=-INFINITY;

  for(ii=0;ii<n;ii++){
    if(abs(tree[0].r[ii*m+0])>tree[0].space[m*0+0]){
      tree[0].space[m*0+0]=abs(tree[0].r[ii*m+0]);
    }
    if(abs(tree[0].r[ii*m+1])>tree[0].space[m*0+1]){
      tree[0].space[m*0+1]=abs(tree[0].r[ii*m+1]);
    }
    if(abs(tree[0].r[ii*m+2])>tree[0].space[m*0+2]){
      tree[0].space[m*0+2]=abs(tree[0].r[ii*m+2]);
    }
  }
  
  /* search for largest dimension */
  if(tree[0].space[m*0+1]<tree[0].space[m*0+2]){
    tree[0].space[m*0+1]=tree[0].space[m*0+2];
  }
  if(tree[0].space[m*0+0]<tree[0].space[m*0+1]){
    tree[0].space[m*0+0]=tree[0].space[m*0+1];
  }

  /* cut some slack */
  tree[0].space[m*0+0]*=2.0;

  /* first corner of root box */
  tree[0].space[m*0+1]=-tree[0].space[m*0+0];
  tree[0].space[m*0+2]=tree[0].space[m*0+0];

  /* second corner of root box */
  tree[0].space[m*1+0]=tree[0].space[m*0+0];
  tree[0].space[m*1+1]=-tree[0].space[m*0+0];
  tree[0].space[m*1+2]=-tree[0].space[m*0+0];
  
  /* third corner of root box */
  tree[0].space[m*2+0]=-tree[0].space[m*0+0];
  tree[0].space[m*2+1]=-tree[0].space[m*0+0];
  tree[0].space[m*2+2]=-tree[0].space[m*0+0];
  
  /* fourth corner of root box */
  tree[0].space[m*3+0]=-tree[0].space[m*0+0];
  tree[0].space[m*3+1]=-tree[0].space[m*0+0];
  tree[0].space[m*3+2]=tree[0].space[m*0+0];
  
  /* fifth corner of root box */
  tree[0].space[m*4+0]=tree[0].space[m*0+0];
  tree[0].space[m*4+1]=tree[0].space[m*0+0];
  tree[0].space[m*4+2]=tree[0].space[m*0+0];
  
  /* sixth corner of root box */
  tree[0].space[m*5+0]=tree[0].space[m*0+0];
  tree[0].space[m*5+1]=tree[0].space[m*0+0];
  tree[0].space[m*5+2]=-tree[0].space[m*0+0];
  
  /* seventh corner of root box */
  tree[0].space[m*6+0]=-tree[0].space[m*0+0];
  tree[0].space[m*6+1]=tree[0].space[m*0+0];
  tree[0].space[m*6+2]=-tree[0].space[m*0+0];
  
  /* eigth corner of root box */
  tree[0].space[m*7+0]=-tree[0].space[m*0+0];
  tree[0].space[m*7+1]=tree[0].space[m*0+0];
  tree[0].space[m*7+2]=tree[0].space[m*0+0];
  
  tree[0].mass=0;
  for(ii=0;ii<n;ii++)
    tree[0].mass+=tree[0].m[ii];
  tree[0].num=n;
  
  tree[0].center[0]=0;
  tree[0].center[1]=0;
  tree[0].center[2]=0;
  for(ii=0;ii<n;ii++){
    tree[0].center[0]+=tree[0].r[ii*m+0]*tree[0].m[ii];
    tree[0].center[1]+=tree[0].r[ii*m+1]*tree[0].m[ii];
    tree[0].center[2]+=tree[0].r[ii*m+2]*tree[0].m[ii];
  }
  tree[0].center[0]/=tree[0].mass;
  tree[0].center[1]/=tree[0].mass;
  tree[0].center[2]/=tree[0].mass;
  tree[0].numcells=1;
  tree[0].parent=0;
  tree[0].index=0;
  tree[0].numchild=0;
  tree[0].l=2*tree[0].space[0*m+0];

  return 0;
}

void direct_summation(struct universe *world, dubfloat_t *r, dubfloat_t *a, dubfloat_t G){
  /* loop variables */
  int ii;

  /* distance calculation variables */
  dubfloat_t d;
  dubfloat_t d2;

  /* plummer softening length */
  dubfloat_t epsilon;

  /* calculate gravitational acceleration by direct summation */
  for(ii=0;ii<world->num;ii++){
    /* compute distance */
    d=euclidean_distance(r, &world->r[ii*3], 3);

    if(d>1E-8){
      epsilon=world->h[ii];
      //d2=sqrt(d*d+epsilon*epsilon);
      //d2=1/(d2*d2*d2);
      d2=kernel_grav_g(d, epsilon);

      a[0]-=(G*world->m[ii])*d2*(r[0]-world->r[ii*3+0]);
      a[1]-=(G*world->m[ii])*d2*(r[1]-world->r[ii*3+1]);
      a[2]-=(G*world->m[ii])*d2*(r[2]-world->r[ii*3+2]);
    }
  }
}

void force_walk(struct universe *world, struct cell *tree, struct cell *root, dubfloat_t *r,
		dubfloat_t *a, dubfloat_t G, dubfloat_t theta, dubfloat_t h){
  /* loop variables */
  int ii;

  /* child cell */
  struct cell *child;

  /* cell node stack */
  int nodestack[tree->numcells];

  /* walk pointer */
  int wp;

  /* cell side length variable */
  dubfloat_t l;

  /* distance calculation variables */
  dubfloat_t d,d2,delta;

  /* plummer softening factor */
  dubfloat_t epsilon;

  /* temporary vector */
  dubfloat_t temp[3];

  /* init walk pointer */
  wp=0;

  /* init stack with root cell siblings */
  for(ii=0;ii<root->numchild;ii++){
    nodestack[wp]=root->children[ii];
    wp++;
  } 

  /* central gravitation */
  /*
  temp[0]=0.0;
  temp[1]=0.0;
  temp[2]=0.0;
  d=euclidean_distance(r, &temp[0], 3);
  d2=kernel_grav_g(d, 1.0);

  a[0]-=(G*0.4)*d2*(r[0]);
  a[1]-=(G*0.4)*d2*(r[1]);
  a[2]-=(G*0.4)*d2*(r[2]);	
  */
  
  while(wp>0){
    /* pop a cell from the node stack */
    wp--;
    child=&tree[nodestack[wp]];

    /* cell side lengths */
    l=child->l;

    /* compute distance */
    d=euclidean_distance(r, child->center, 3);

    /* compute distance from cells geometrical center to center of mass */
    /* for cell opening criterion */
    temp[0]=child->center[0]-(child->space[0*3+0]-(l/2.0));
    temp[1]=child->center[1]-(child->space[0*3+1]+(l/2.0));
    temp[2]=child->center[2]-(child->space[0*3+2]-(l/2.0));
    delta=euclidean_norm(&temp[0], 3);

    if(child->num>1){
      if(d>(child->l/theta)+delta){
	/* approximate cell as ensemble */
	/* calculate acceleration with softened potential */
	epsilon=child->distr_len;	  
	d2=kernel_grav_g(d, epsilon);

	a[0]-=(G*child->mass)*d2*(r[0]-child->center[0]);
	a[1]-=(G*child->mass)*d2*(r[1]-child->center[1]);
	a[2]-=(G*child->mass)*d2*(r[2]-child->center[2]);	
      }
      else{
	/* push sibling nodes on the stack */
	for(ii=0;ii<child->numchild;ii++){
	  nodestack[wp]=child->children[ii];
	  wp++;
	}
      }
    }
    else{
      /* single particle cell */
      /* calculate acceleration with softened potential */
      epsilon=world->h[child->particle_index];
      d2=kernel_grav_g(d, epsilon);

      a[0]-=(G*child->mass)*d2*(r[0]-child->center[0]);
      a[1]-=(G*child->mass)*d2*(r[1]-child->center[1]);
      a[2]-=(G*child->mass)*d2*(r[2]-child->center[2]);
    }
  }
}

void neighbour_walk(struct cell *tree, struct cell *root, dubfloat_t *r, dubfloat_t h, dubfloat_t max_h,
		    dubfloat_t *h_in, int *neighbour_num, int *neighbour_list){
  /* loop variables */
  int ii;

  /* child cell */
  struct cell *child;

  /* cell node stack */
  int nodestack[tree->numcells];

  /* walk pointer */
  int wp;

  /* cell side length variable */
  dubfloat_t l;

  /* distance calculation variables */
  dubfloat_t d;

  /* temporary vector */
  dubfloat_t temp[3];

  /* init walk pointer */
  wp=0;

  /* init stack with root cell siblings */
  for(ii=0;ii<root->numchild;ii++){
    nodestack[wp]=root->children[ii];
    wp++;
  } 

  while(wp>0){
    /* pop a cell from the node stack */
    wp--;
    child=&tree[nodestack[wp]];

    /* cell side lengths */
    l=child->l;

    if(child->num>1){
      /* compute distance from cells geometrical center to particle */
      temp[0]=r[0]-(child->space[0*3+0]-l/2.0);
      temp[1]=r[1]-(child->space[0*3+1]+l/2.0);
      temp[2]=r[2]-(child->space[0*3+2]-l/2.0);
      d=euclidean_norm(&temp[0],3);

      /* if the smoothing length radius intersects with the cells corner */
      /* radius, recurse deeper into the tree to find more particles */
      if((sqrt(3.0)/2.0)*l+max_h+2.0*h>d){
	/* push sibling nodes into the stack */
	for(ii=0;ii<child->numchild;ii++){
	  nodestack[wp]=child->children[ii];
	  wp++;
	}
      }
    }
    else{
      /* compute distance from cell center of mass to particle */
      d=euclidean_distance(&r[0],&child->center[0],3);

      if(d/h<2.0||d/h_in[child->particle_index]<2.0){
	/* add particle index into the neighbouring particle list */
	neighbour_list[*neighbour_num]=child->particle_index;
	*neighbour_num+=1;
      }
    }
  }
}

void neighbour_recurse(struct cell *tree, struct cell *root, dubfloat_t *r, dubfloat_t h, dubfloat_t max_h,
		       dubfloat_t *h_in, int *neighbour_num, int *neighbour_list){
  /* loop variables */
  int ii;

  /* distance calculation variable */
  dubfloat_t d;

  /* cell length variable */
  dubfloat_t l;

  /* temporary vector */
  dubfloat_t temp[3];

  /* tree structure child pointer */
  struct cell *child;
    
  /* check thru child cells */
  for(ii=0;ii<root->numchild;ii++){
    /* child cell */
    child=&tree[root->children[ii]];

    /* cell side lengths */
    l=child->l;

    if(child->num>1){
      /* compute distance from cell centers to particles */
      temp[0]=r[0]-(child->space[0*3+0]-l/2.0);
      temp[1]=r[1]-(child->space[0*3+1]+l/2.0);
      temp[2]=r[2]-(child->space[0*3+2]-l/2.0);
      d=euclidean_norm(&temp[0],3);

      /* if the smoothing length radius intersects with the cells corner */
      /* radius, recurse deeper into the tree to find more particles */
      if((sqrt(3.0)/2.0)*l+max_h+2.0*h>d){
	neighbour_recurse(tree, child, r, h, max_h, h_in, neighbour_num, neighbour_list);
      }
    }
    else{
      /* compute distance from cell center of mass to particle */
      d=euclidean_distance(&r[0],&child->center[0],3);

      if(d/h<2.0||d/h_in[child->particle_index]<2.0){
	/* add particle index into the neighbouring particle list */
	neighbour_list[*neighbour_num]=child->particle_index;
	*neighbour_num+=1;
      }
    }
  }
}

void compute_total_energy(struct universe *world, dubfloat_t theta, int lo, int hi){
  /* loop variables */
  int ii;

  /* pointers to state vectors */
  dubfloat_t *r_in;
  dubfloat_t *m_in;

  /* total gravitational potential energy */
  dubfloat_t *U;

  /* total kinetic energy */
  dubfloat_t *K;

  /* pointer to barnes hut tree */
  struct cell *tree;

  r_in=world->r;
  m_in=world->m;

  tree=world->tree;

  U=&world->u_grav;
  K=&world->u_kin;

  /* sum total potential and kinetic energy */
  for(ii=lo;ii<hi;ii++){
    /* sum in potential energy from tree */
    /* unit is AU^3/(M_solar*yr^2)*M_solar^2/AU */
    potential_recurse(tree, &tree[0], &r_in[3*ii], m_in[ii], U, world->G, theta, world->h[ii]);
    //potential_recurse(tree, &tree[0], &r_in[3*ii], m_in[ii], U, world->G, theta, world->epsilon);

    /* sum in kinetic energy */
    /* unit is M_solar*(AU/yr)^2 */
    *K+=0.5*world->m[ii]*(world->v[3*ii+0]*world->v[3*ii+0]+world->v[3*ii+1]*world->v[3*ii+1]+
			  world->v[3*ii+2]*world->v[3*ii+2]);
  }
}

void potential_recurse(struct cell *tree, struct cell *root, dubfloat_t *r, dubfloat_t m,
		       dubfloat_t *U, dubfloat_t G, dubfloat_t theta, dubfloat_t epsilon){
  /* loop variables */
  int ii;

  /* distance calculation variables */
  dubfloat_t d;
  dubfloat_t d2;
  //dubfloat_t delta;

  /* cell length variable */
  //dubfloat_t l;

  /* temporary vector */
  //dubfloat_t temp[3];

  /* tree structure child pointer */
  struct cell *child;
    
  /* check thru child cells */
  for(ii=0;ii<root->numchild;ii++){
    child=&tree[root->children[ii]];

    /* cell side lengths */
    //l=child->l;

    /* compute distance */
    d=euclidean_distance(&r[0],&child->center[0],3);
    
    /* compute distance from cell geometrical center to center of mass */
    /*
    temp[0]=child->center[0]-(child->space[0*3+0]-l/2.0);
    temp[1]=child->center[1]-(child->space[0*3+1]+l/2.0);
    temp[2]=child->center[2]-(child->space[0*3+2]-l/2.0);
    delta=euclidean_norm(&temp[0],3);
    */

    if(child->num>1){
      if(0){
      //if(d>child->l/theta+delta){
	/* approximate as ensemble */
	/* calculate potential energy with plummer softening */
	epsilon=child->distr_len;	  
	//d2=1/sqrt(d*d+epsilon*epsilon);	
	d2=kernel_grav_f(d, epsilon);
	*U-=G*m*child->mass*d2;
      }
      else{
	/* recurse deeper into tree */
	potential_recurse(tree, child, r, m, U, G, theta, epsilon);
      }
    }
    else{
      /* calculate potential energy with plummer softening */
      epsilon=child->distr_len;	  
      d2=kernel_grav_f(d, epsilon);
      *U-=G*m*child->mass*d2;
    }
  }
}

void force_recurse(struct cell *tree, struct cell *root, dubfloat_t *r,
		   dubfloat_t *f, dubfloat_t G, dubfloat_t theta, dubfloat_t epsilon){
  /* loop variables */
  int ii;

  /* distance calculation variables */
  dubfloat_t d;
  dubfloat_t d2;
  dubfloat_t delta;

  /* cell length variable */
  dubfloat_t l;

  /* temporary vector */
  dubfloat_t temp[3];

  /* tree structure child pointer */
  struct cell *child;
    
  /* check thru child cells */
  for(ii=0;ii<root->numchild;ii++){
    child=&tree[root->children[ii]];

    /* cell side lengths */
    l=child->l;

    /* compute distance */
    d=euclidean_distance(&r[0],&child->center[0],3);
    
    /* compute distance from cell geometrical center to center of mass */
    temp[0]=child->center[0]-(child->space[0*3+0]-l/2.0);
    temp[1]=child->center[1]-(child->space[0*3+1]+l/2.0);
    temp[2]=child->center[2]-(child->space[0*3+2]-l/2.0);
    delta=euclidean_norm(&temp[0],3);

    if(child->num>1){
      if(d>child->l/theta+delta){
	/* approximate as ensemble */
	/* calculate force with softened potential */
	d2=kernel_grav_g(d, epsilon);
	f[0]-=(G*child->mass)*d2*(r[0]-child->center[0]);
	f[1]-=(G*child->mass)*d2*(r[1]-child->center[1]);
	f[2]-=(G*child->mass)*d2*(r[2]-child->center[2]);	
      }
      else{
	/* recurse deeper into tree */
	force_recurse(tree, child, r, f, G, theta, epsilon);
      }
    }
    else{
      /* calculate force with softened potential */
      d2=kernel_grav_g(d, epsilon);
      f[0]-=(G*child->mass)*d2*(r[0]-child->center[0]);
      f[1]-=(G*child->mass)*d2*(r[1]-child->center[1]);
      f[2]-=(G*child->mass)*d2*(r[2]-child->center[2]);
    }
  }
}

void tree_recurse(struct cell *tree, struct cell *root){
    int ii;

    printf("\nroot->index: %d\n", root->index);
    printf("root->num: %d\n", root->num);
    printf("root->mass: %f\n", root->mass);    
    printf("root->center[0]: %f\n", root->center[0]);
    printf("root->center[1]: %f\n", root->center[1]);
    printf("root->center[2]: %f\n", root->center[2]);
    
    printf("root->space: \n");
    for(ii=0;ii<8;ii++)
      printf("%f %f %f\n", root->space[ii*3+0], root->space[ii*3+1], root->space[ii*3+1]);

    printf("root->numchild: %d\n", root->numchild);

    printf("root->children: ");
    for(ii=0;ii<root->numchild;ii++)
        printf("%d ", root->children[ii]);
    printf("\n");

    printf("root->r: ");
    for(ii=0;ii<root->num;ii++)
      printf("%f %f %f\n", root->r[ii*3+0], root->r[ii*3+1], root->r[ii*3+2]);
    printf("\n");

    if(root->num>1){
        for(ii=0;ii<root->numchild;ii++){
            tree_recurse(tree, &tree[root->children[ii]]);
        }
    }
}

void branch_recurse(struct universe *world, struct cell *tree, struct cell *root, int *cellindex){
  int ii;
  struct cell *child;

  /* branch into subcells if root cell has more than one particle */
  if(root->num>1){
    tree_branch(world, tree, root, cellindex);
    free(root->r);
    free(root->m);
    free(root->particle_index_list);

    for(ii=0;ii<root->numchild;ii++){
      child=&tree[root->children[ii]];
      branch_recurse(world, tree, child, cellindex);
    }
  }
}

void tree_branch(struct universe *world, struct cell *tree, struct cell *root, int *cellindex){
    struct cell *newcell;
    dubfloat_t x,y,z;
    dubfloat_t *r,*m;
    int *particle_index_list;
    dubfloat_t h;
    dubfloat_t cell_r[3*8];
    int ii,jj,nn;
    dubfloat_t hh;

    /* root particle displacement, mass and particle index vectors */
    r=root->r;
    m=root->m;
    particle_index_list=root->particle_index_list;
    
    /* root cell side length divided by two to get the subcell side length */
    h=(root->space[0*3+2]-root->space[1*3+2])/2.0;

    /* subcell corner displacement vectors */
    cell_r[0*3+0]=root->space[0*3+0];
    cell_r[0*3+1]=root->space[0*3+1];
    cell_r[0*3+2]=root->space[0*3+2];

    cell_r[1*3+0]=root->space[0*3+0];
    cell_r[1*3+1]=root->space[0*3+1];
    cell_r[1*3+2]=root->space[0*3+2]-h;

    cell_r[2*3+0]=root->space[0*3+0]-h;
    cell_r[2*3+1]=root->space[0*3+1];
    cell_r[2*3+2]=root->space[0*3+2]-h;

    cell_r[3*3+0]=root->space[0*3+0]-h;
    cell_r[3*3+1]=root->space[0*3+1];
    cell_r[3*3+2]=root->space[0*3+2];

    cell_r[4*3+0]=root->space[0*3+0];
    cell_r[4*3+1]=root->space[0*3+1]+h;
    cell_r[4*3+2]=root->space[0*3+2];

    cell_r[5*3+0]=root->space[0*3+0];
    cell_r[5*3+1]=root->space[0*3+1]+h;
    cell_r[5*3+2]=root->space[0*3+2]-h;

    cell_r[6*3+0]=root->space[0*3+0]-h;
    cell_r[6*3+1]=root->space[0*3+1]+h;
    cell_r[6*3+2]=root->space[0*3+2]-h;

    cell_r[7*3+0]=root->space[0*3+0]-h;
    cell_r[7*3+1]=root->space[0*3+1]+h;
    cell_r[7*3+2]=root->space[0*3+2];

    /*
    cell_r[0*2+0]=root->space[0*2+0];
    cell_r[0*2+1]=root->space[0*2+1];
    cell_r[1*2+0]=root->space[0*2+0];
    cell_r[1*2+1]=root->space[0*2+1]-h;
    cell_r[2*2+0]=root->space[0*2+0]-h;
    cell_r[2*2+1]=root->space[0*2+1]-h;
    cell_r[3*2+0]=root->space[0*2+0]-h;
    cell_r[3*2+1]=root->space[0*2+1];
    */

    /* form cell subcells*/
    for(jj=0;jj<8;jj++){
      /* add new subcell candidate to the tree */
      newcell=&tree[tree[0].numcells];
      newcell->index=tree[0].numcells;

      /* subcell corner displacement vector */
      x=cell_r[jj*3+0];
      y=cell_r[jj*3+1];
      z=cell_r[jj*3+2];
    
      /* form subcells displacement vectors from corner displacement */
      newcell->space[0*3+0]=x;
      newcell->space[0*3+1]=y;
      newcell->space[0*3+2]=z;

      newcell->space[1*3+0]=x;
      newcell->space[1*3+1]=y;
      newcell->space[1*3+2]=z-h;

      newcell->space[2*3+0]=x-h;
      newcell->space[2*3+1]=y;
      newcell->space[2*3+2]=z-h;

      newcell->space[3*3+0]=x-h;
      newcell->space[3*3+1]=y;
      newcell->space[3*3+2]=z;

      newcell->space[4*3+0]=x;
      newcell->space[4*3+1]=y+h;
      newcell->space[4*3+2]=z;

      newcell->space[5*3+0]=x;
      newcell->space[5*3+1]=y+h;
      newcell->space[5*3+2]=z-h;

      newcell->space[6*3+0]=x-h;
      newcell->space[6*3+1]=y+h;
      newcell->space[6*3+2]=z-h;

      newcell->space[7*3+0]=x-h;
      newcell->space[7*3+1]=y+h;
      newcell->space[7*3+2]=z;

      /* reset subcell variables */
      newcell->center[0]=0;
      newcell->center[1]=0;
      newcell->center[2]=0;
      newcell->mass=0;
      newcell->num=0;
      newcell->numchild=0;
      newcell->l=h;

      /* allocate memory for particle displacement, mass and index vectors */
      newcell->r=(dubfloat_t*)malloc(3*root->num*sizeof(dubfloat_t));
      newcell->m=(dubfloat_t*)malloc(root->num*sizeof(dubfloat_t));
      newcell->particle_index_list=(int*)malloc(root->num*sizeof(int));

      /* calculate new subcell particle displacement and mass vectors */
      nn=0;
      for(ii=0;ii<root->num;ii++){
        if(r[ii*3+0]<newcell->space[0*3+0]&&r[ii*3+1]<newcell->space[4*3+1]&&
           r[ii*3+0]>newcell->space[3*3+0]&&r[ii*3+1]>newcell->space[0*3+1]&&
           r[ii*3+2]<newcell->space[0*3+2]&&r[ii*3+2]>newcell->space[1*3+2]
	  ){
	  newcell->r[nn*3+0]=r[ii*3+0];
	  newcell->r[nn*3+1]=r[ii*3+1];
	  newcell->r[nn*3+2]=r[ii*3+2];
	  newcell->m[nn]=m[ii];
	  newcell->particle_index_list[nn]=particle_index_list[ii];
	  newcell->center[0]+=m[ii]*r[ii*3+0];
	  newcell->center[1]+=m[ii]*r[ii*3+1];
	  newcell->center[2]+=m[ii]*r[ii*3+2];
	  newcell->mass+=m[ii];
	  newcell->num++;
	  nn++;
        }
      }
      /* subcell parent cell */
      newcell->parent=root->index;

      /* add valid subcell to the tree */
      if(newcell->num>0){
        newcell->center[0]/=newcell->mass;
        newcell->center[1]/=newcell->mass;
        newcell->center[2]/=newcell->mass;
        tree[0].numcells++;
        root->children[root->numchild++]=newcell->index;

	/* approximate mass distribution */
	/* by averaging gravitational softening lengths */
	hh=0;
	for(ii=0;ii<newcell->num;ii++){
	  hh+=world->h[newcell->particle_index_list[ii]];
	}
	hh/=newcell->num;
	newcell->distr_len=hh;

	/* write particle cell indeces if we reached the end of a branch */
	if(newcell->num==1){
	  newcell->particle_index=newcell->particle_index_list[0];
	  cellindex[newcell->particle_index]=newcell->index;
	}
	else{
	  newcell->particle_index=-1;
	}
      }

      /* if subcell has only 1 or no particles, free useless vectors */
      /* useless as subcell wont have child cells */
      if(newcell->num<2){
        free(newcell->r);
        free(newcell->m);
        free(newcell->particle_index_list);
      }
    }
}
