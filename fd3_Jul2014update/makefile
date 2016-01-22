
all: fd3 examples archive

clean :
	rm fd3 *.o art_??????.??? V453_Cyg.??? fd3.tar.gz

fd3 : fd3fpolis.o fd3sep.o triorb.o kepler.o mxfuns.o
	${CC} -Wall fd3fpolis.o fd3sep.o triorb.o kepler.o mxfuns.o \
		-lgsl -lgslcblas -lm -o $@

.c.o :
	${CC} -Wall -c $<

examples :
	./fd3 < art_single.in > art_single.out
	./fd3 < art_double.in > art_double.out
	./fd3 < art_triple.in > art_triple.out
	./fd3 < V453_Cyg.in > V453_Cyg.out

archive :
	tar -cvf fd3.tar *
	gzip fd3.tar

