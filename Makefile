# compile flags
THREAD_PROFILING=0
ENERGY_PROFILING=0
CUDA=0
SDL=1
OPTIMIZATION_LEVEL=3

# object files
ifeq ($(CUDA), 1)
	OBJ=dubstep.o tree.o sph.o threads.o timer.o ic.o lattice.o statistics.o sph_cuda.o
else
	OBJ=dubstep.o tree.o sph.o threads.o timer.o ic.o lattice.o statistics.o 
endif

# compilers
CC=gcc
NVCC=/usr/local/cuda-5.0/bin/nvcc

# compiler options
ifeq ($(SDL),1)
	OPTS=-Wall -pthread `sdl-config --cflags`
else
	OPTS=-Wall -pthread
endif

# compiler flags
ifeq ($(CUDA), 1)
	ifeq ($(SDL),1)
		CFLAGS=-DCUDA=1 -DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -finline-functions -I/usr/local/cuda-5.0/include -lm -lGL -lrt -ltiff `sdl-config --libs` -L/usr/local/cuda-5.0/lib -lcudart
	else
		CFLAGS=-DCUDA=1 -DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -finline-functions -I/usr/local/cuda-5.0/include -lm -lGL -lrt -L/usr/local/cuda-5.0/lib -lcudart
	endif
else
	ifeq ($(SDL),1)
		CFLAGS=-DCUDA=0 -DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -lm -lGL -ltiff -lrt `sdl-config --libs`
	else
		CFLAGS=-DCUDA=0 -DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -lm -lrt
	endif
endif

ifeq ($(CUDA), 1)
dubstep: $(OBJ)
	$(NVCC) -m32 -o $@ $+ -DCUDA=1 -DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -I/usr/local/cuda-5.0/include -lm -lGL -lrt -ltiff `sdl-config --libs` -L/usr/local/cuda-5.0/lib -lcudart
else
dubstep: $(OBJ)
	$(CC) -o $@ $+ $(OPTS) $(CFLAGS) 
endif

dubstep.o: dubstep.c
	$(CC) $(OPTS) $(CFLAGS) -c $<

tree.o: tree.c tree.h
	$(CC) $(OPTS) $(CFLAGS) -c $<

sph.o: 	sph.c sph.h
	$(CC) $(OPTS) $(CFLAGS) -c $<

ifeq ($(CUDA), 1)
sph_cuda.o: sph_cuda.cu sph_cuda.h
	$(NVCC) -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -m32 -I. -I/usr/local/cuda-5.0/include -c $<
endif

threads.o: threads.c threads.h
	$(CC) $(OPTS) $(CFLAGS) -c $<

timer.o: timer.c timer.h
	$(CC) $(OPTS) $(CFLAGS) -c $<

ic.o: ic.c ic.h
	$(CC) $(OPTS) $(CFLAGS) -c $<

lattice.o: lattice.c lattice.h
	$(CC) $(OPTS) $(CFLAGS) -c $<

statistics.o: statistics.c statistics.h
	$(CC) $(OPTS) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm *.o dubstep

