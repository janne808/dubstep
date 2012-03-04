#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dubstep.h"
#include "tree.h"

int init_treeroot(struct cell *tree, struct universe *world){
  int ii;
  int n,m;
  
  /* determine state vector dimensions */
  n=world->num;
  m=world->dim;

  /* form tree root cell*/
  /* init particle displacement vector */
  tree[0].r=malloc(m*n*sizeof(double));
  if(!tree[0].r){
    printf("Out of memory: root displacement vector not allocated.\n");
    return -1;
  }
  memcpy(tree[0].r, world->r , m*n*sizeof(double));

  /* init particle mass vector */
  tree[0].m=malloc(n*sizeof(double));        
  if(!tree[0].m){
    printf("Out of memory: root mass vector not allocated.\n");
    return -1;
  }
  memcpy(tree[0].m, world->m , n*sizeof(double));
  
  /* init particle index vector */
  /*
  tree[0].particle_index=malloc(n*sizeof(int));
  if(!tree[0].particle_index){
    printf("Out of memory: root particle index vector not allocated.\n");
    return -1;
  }
  for(ii=0;ii<n;ii++){
    tree[0].particle_index[ii]=ii;
  }
  */

  /* search for maximum values */
  /* make sure we are dealing with squares */
  tree[0].space[m*0+0]=-INFINITY;
  tree[0].space[m*0+1]=-INFINITY;
  tree[0].space[m*0+2]=-INFINITY;
  for(ii=0;ii<n;ii++){
    if(abs(tree[0].r[ii*m+0])>tree[0].space[m*0+0]){
      tree[0].space[m*0+0]=tree[0].r[ii*m+0];
    }
    if(abs(tree[0].r[ii*m+1])>tree[0].space[m*0+1]){
      tree[0].space[m*0+1]=tree[0].r[ii*m+1];
    }
    if(abs(tree[0].r[ii*m+2])>tree[0].space[m*0+2]){
      tree[0].space[m*0+2]=tree[0].r[ii*m+2];
    }
  }
  
  /* search for largest coordinate */
  if(tree[0].space[m*0+1]<tree[0].space[m*0+2]){
    tree[0].space[m*0+1]=tree[0].space[m*0+2];
  }
  if(tree[0].space[m*0+0]<tree[0].space[m*0+1]){
    tree[0].space[m*0+0]=tree[0].space[m*0+1];
  }

  tree[0].space[m*0+0]+=0.5;

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

void forcerecurse(struct cell *tree, struct cell *root, double *r,
		  double *f, double G, double theta, double epsilon){
  int ii;
  double dx,dy,dz,d,d2;
  struct cell *child;
    
  /* check thru child cells */
  for(ii=0;ii<root->numchild;ii++){
    child=&tree[root->children[ii]];

    /* compute distance */
    dx=r[0]-child->center[0];
    dy=r[1]-child->center[1];
    dz=r[2]-child->center[2];
    d=sqrt(dx*dx+dy*dy+dz*dz);
    
    if(child->num>1){
      if(child->l/d>theta){
	/* open new cell */
	forcerecurse(tree, child, r, f, G, theta, epsilon);
      }
      else{
	/* approximate as ensemble */
	/* calculate force with plummer softening */
	d2=sqrt(d*d+epsilon*epsilon);
	d2=1/(d2*d2*d2);
	f[0]-=(G*child->mass)*d2*(r[0]-child->center[0]);
	f[1]-=(G*child->mass)*d2*(r[1]-child->center[1]);
	f[2]-=(G*child->mass)*d2*(r[2]-child->center[2]);	
      }
    }
    else{
      /* calculate force with plummer softening */
      d2=sqrt(d*d+epsilon*epsilon);
      d2=1/(d2*d2*d2);
      f[0]-=(G*child->mass)*d2*(r[0]-child->center[0]);
      f[1]-=(G*child->mass)*d2*(r[1]-child->center[1]);
      f[2]-=(G*child->mass)*d2*(r[2]-child->center[2]);
    }
  }
}

void treerecurse(struct cell *tree, struct cell *root){
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
            treerecurse(tree, &tree[root->children[ii]]);
        }
    }
}

void branchrecurse(struct cell *tree, struct cell *root, int *cellindex){
  int ii;
  if(root->num>1){
    treebranch(tree, root, cellindex);
    free(root->r);
    free(root->m);
    //free(root->particle_index);
    for(ii=0;ii<root->numchild;ii++)
      branchrecurse(tree, &tree[root->children[ii]], cellindex);
  }
}

void treebranch(struct cell *tree, struct cell *root, int *cellindex){
    struct cell *newcell;
    double x,y,z;
    double *r,*m;
    //int *particle_index;
    double h;
    double cell_r[3*8];
    int ii,jj,nn;

    /* root particle displacement, mass and particle index vectors */
    r=root->r;
    m=root->m;
    //particle_index=root->particle_index;
    
    /* root cell side length divided by 2 */
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
      newcell->r=(double*)malloc(3*root->num*sizeof(double));
      newcell->m=(double*)malloc(root->num*sizeof(double));
      //newcell->particle_index=(int*)malloc(root->num*sizeof(int));

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
	  //newcell->particle_index[nn]=particle_index[ii];
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

	/* write particle cell index if we reached the end of a branch */
	//if(newcell->num==1)
	//  cellindex[newcell->particle_index[0]]=newcell->index;
      }

      /* if subcell has only 1 or no particles, free useless vectors */
      /* useless as subcell wont have child cells*/
      if(newcell->num<2){
        free(newcell->r);
        free(newcell->m);
      }
    }
}