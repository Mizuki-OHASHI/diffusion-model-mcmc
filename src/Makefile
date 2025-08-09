CC = gcc
# CFLAGS = -O3
CFLAGS = -O3 -fopenmp -Wall -Wextra -Wpedantic -std=c99

default: ising2d_mh_run ising2d_mh_sweep xy2d_mh_run xy2d_mh_sweep

utils/stat.o: utils/stat.c utils/stat.h
	$(CC) $(CFLAGS) -c utils/stat.c -o utils/stat.o

ising2d_mh.o: ising2d_mh.c ising2d_mh.h utils/stat.h
	$(CC) $(CFLAGS) -c ising2d_mh.c

ising2d_mh_run: ising2d_mh_run.c ising2d_mh.o utils/stat.o
	$(CC) $(CFLAGS) -o ising2d_mh_run.out ising2d_mh_run.c ising2d_mh.o utils/stat.o -lm


ising2d_mh_sweep: ising2d_mh_sweep.c ising2d_mh.o utils/stat.o
	$(CC) $(CFLAGS) -o ising2d_mh_sweep.out ising2d_mh_sweep.c ising2d_mh.o utils/stat.o -lm

xy2d_mh.o: xy2d_mh.c xy2d_mh.h utils/stat.h
	$(CC) $(CFLAGS) -c xy2d_mh.c

xy2d_mh_run: xy2d_mh_run.c xy2d_mh.o utils/stat.o
	$(CC) $(CFLAGS) -o xy2d_mh_run.out xy2d_mh_run.c xy2d_mh.o utils/stat.o -lm

xy2d_mh_sweep: xy2d_mh_sweep.c xy2d_mh.o utils/stat.o
	$(CC) $(CFLAGS) -o xy2d_mh_sweep.out xy2d_mh_sweep.c xy2d_mh.o utils/stat.o -lm

mv:
	mv *.dat ../data/

clean:
	rm -f *.o *.out
