
OBJ=dubstep.o tree.o sph.o threads.o
CC=gcc
OPTS=-Wall -g -pthread `sdl-config --cflags`
CFLAGS=-DINFINITY=HUGE_VAL -DENABLE_GUI=1 -finline-functions -lm -lGL -ltiff `sdl-config --libs`

dubstep: $(OBJ)
	$(CC) -o $@ $(OPTS) $(CFLAGS) $(OBJ) 

dubstep.o: dubstep.c tree.h
	gcc $(OPTS) $(CFLAGS) -c $<

tree.o: tree.c tree.h
	gcc $(OPTS) $(CFLAGS) -c $<

sph.o: 	sph.c sph.h
	gcc $(OPTS) $(CFLAGS) -c $<

threads.o: threads.c threads.h
	gcc $(OPTS) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm *.o dubstep

