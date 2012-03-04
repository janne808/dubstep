#ifndef __TREEH__
#define __TREEH__

/* tree cell */
struct cell{ 
  double *r; /* cell particle displacement vector */
  double *m; /* cell particle mass vector */

  int *particle_index; /* cell particle index relative to the world vector */

  double space[3*8]; /* cell spatial limits */
    
  double l;   /* cell size */
    
  double mass; /* total mass of particles */
  int num; /* number of particles in cell */
  double center[3]; /* cell center of mass */
    
  int parent; /* parent cell index */
  int children[8]; /* children cell indices */
  int numcells; /* amount of cells in tree */
  int index;
  int numchild;
};

/* function prototypes */
int init_treeroot(struct cell *tree, struct universe *world);
void treebranch(struct cell *tree, struct cell *root, int *cellindex);
void treerecurse(struct cell *tree, struct cell *root);
void forcerecurse(struct cell *tree, struct cell *root,
		  double *r, double *f, double G, double theta,
		  double epsilon);
void branchrecurse(struct cell *tree, struct cell *root, int *cellindex);

#endif