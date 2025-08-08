CC = gcc
CFLAGS = -O3

default: ising2d_mh

ising2d_mh: ising2d_mh.c utils/stat.c
	$(CC) $(CFLAGS) -o ising2d_mh.out ising2d_mh.c utils/stat.c -lm

mv:
	mv ising2d_mh*.dat ../data/

clean:
	rm -f ising2d_mh.out
