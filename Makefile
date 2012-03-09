
OBJ=dubstep.o tree.o sph.o threads.o sph_cuda.o
CC=gcc
NVCC=nvcc
OPTS=-Wall -g -pthread `sdl-config --cflags`
CFLAGS=-DINFINITY=HUGE_VAL -DENABLE_GUI=0 -finline-functions -lm -lGL -ltiff `sdl-config --libs` -L/opt/cuda/lib -lcudart -L/opt/cuda/sdk/lib -lcutil

dubstep: $(OBJ)
	$(CC) -o $@ $(OPTS) $(CFLAGS) $(OBJ) 

dubstep.o: dubstep.c tree.h
	gcc $(OPTS) $(CFLAGS) -c $<

tree.o: tree.c tree.h
	gcc $(OPTS) $(CFLAGS) -c $<

sph.o: 	sph.c sph.h
	gcc $(OPTS) $(CFLAGS) -c $<

sph_cuda.o: sph_cuda.cu sph_cuda.h
	nvcc --compiler-options -fno-strict-aliasing -I. -I/opt/cuda/NVIDIA_CUDA_SDK/common/inc -I/opt/cuda/include -c $<

threads.o: threads.c threads.h
	gcc $(OPTS) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm *.o dubstep

