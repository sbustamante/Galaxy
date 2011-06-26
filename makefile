CC=gcc
CFLAGS=-g -I. -c

#Principal Program
Galaxy.out:data.o numeric.o initial.o galaxy.o potentials.o
	gcc -lm data.o numeric.o initial.o galaxy.o potentials.o -o Galaxy.out
	rm -r *.o

edit:
	kate *.c *.h &

clean:
	rm -r *.o *.out *.png *.tmp script.gpl

cleandatas:
	rm ./datas/*