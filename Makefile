# compile flags
THREAD_PROFILING=0
ENERGY_PROFILING=0
CUDA=0
SDL=0
OPTIMIZATION_LEVEL=3

# object files
ifeq ($(CUDA), 1)
	OBJ=dubstep.o tree.o sph.o threads.o timer.o sph_cuda.o
else
	OBJ=dubstep.o tree.o sph.o threads.o timer.o
endif

# compilers
CC=gcc
NVCC=nvcc

# compiler options
OPTS=-Wall -pthread `sdl-config --cflags`

# compiler flags
ifeq ($(CUDA), 1)
	ifeq ($(SDL),1)
		CFLAGS=-DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -finline-functions -lm -lGL -ltiff `sdl-config --libs` -L/opt/cuda/lib -lcudart -L/opt/cuda/sdk/lib -lcutil
	else
		CFLAGS=-DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -finline-functions -lm -lGL -ltiff -L/opt/cuda/lib -lcudart -L/opt/cuda/sdk/lib -lcutil
	endif
else
	ifeq ($(SDL),1)
		CFLAGS=-DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -lm -lGL -ltiff -lrt `sdl-config --libs`
	else
		CFLAGS=-DINFINITY=HUGE_VAL -DENABLE_GUI=$(SDL) -DTHREAD_PROFILING=$(THREAD_PROFILING) -DENERGY_PROFILING=$(ENERGY_PROFILING) -O$(OPTIMIZATION_LEVEL) -lm -lGL -ltiff -lrt
	endif
endif

dubstep: $(OBJ)
	$(CC) -o $@ $(OPTS) $(CFLAGS) $(OBJ) 

dubstep.o: dubstep.c tree.h
	gcc $(OPTS) $(CFLAGS) -c $<

tree.o: tree.c tree.h
	gcc $(OPTS) $(CFLAGS) -c $<

sph.o: 	sph.c sph.h
	gcc $(OPTS) $(CFLAGS) -c $<

ifeq ($(CUDA), 1)
sph_cuda.o: sph_cuda.cu sph_cuda.h
	nvcc --compiler-options -fno-strict-aliasing -I. -I/opt/cuda/NVIDIA_CUDA_SDK/common/inc -I/opt/cuda/include -c $<
endif

threads.o: threads.c threads.h
	gcc $(OPTS) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm *.o dubstep

