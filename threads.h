#ifndef __THREADSH__
#define __THREADSH__

struct thread_data{
  int thread_id;
  struct universe *world;
  double var;
  int lo;
  int hi;
};

struct thread_data2{
  int thread_id;
  struct universe *world;
  double *r;
  double *v;
  double *a;
  int lo;
  int hi;
};

struct thread_data3{
  int thread_id;
  struct universe *world;
  struct cell *tree;
  double *a_sph;
  double var1;
  double var2;
  double var3;
  int lo;
  int hi;
};

struct thread_data4{
  int thread_id;
  struct universe *world;
  int var1;
  int var2;
  int lo;
  int hi;
};

struct thread_data5{
  int thread_id;
  struct universe *world;
  struct cell *tree;
  double *a_sph;
  double var1;
  double var2;
  double var3;
  int lo;
  int hi;
};

void *smoothingThread(void *threadarg);
void *densityThread(void *threadarg);
void *pressureThread(void *threadarg);
void *soundspeedThread(void *threadarg);
void *CFLThread(void *threadarg);
void *energyThread(void *threadarg);
void *accelerationThread(void *threadarg);
void *integrationThread(void *threadarg);
void *predictorThread(void *threadarg);
void *correctorThread(void *threadarg);

void createSmoothingThreads(struct universe *world, int iterations, int neighbours);
void createDensityThreads(struct universe *world);
void createPressureThreads(struct universe *world);
void createSoundspeedThreads(struct universe *world);
void createCFLThreads(struct universe *world);
void createEnergyThreads(struct universe *world);
void createAccelerationThreads(struct universe *world);
void createIntegrationThreads(struct universe *world);
void createPredictorThreads(struct universe *world);
void createCorrectorThreads(struct universe *world);

#endif
