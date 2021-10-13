CC      = icpc
FC	= ifort

CFLAGS	= -I/u/hikage/work/software/fftw-2.1.5/rfftw -g -I/u/hikage/work/software/fftw-2.1.5/fftw -g

LFLAGS	= -L/u/hikage/work/software/fftw-2.1.5/rfftw/.libs -L/u/hikage/work/software/fftw-2.1.5/fftw/.libs

LIBS	= -lrfftw -lfftw

PROG	= ellipfit

SUB	= main.o fftlog.o

.SUFFIXES: .cpp

.cpp.o	:
	$(CC) -c $(CFLAGS) $<

default	: $(PROG)

ellipfit: $(SUB)
	$(CC) $^ $(LFLAGS) $(LIBS) -o $@

fftlog.o: fftlog.f
	$(FC) -c $^ -nofor_main

test	:
	./ellipfit -i1 -e0.3 -f30 -cimg.dat

random  :
	./ellipfit -i1 -ay -n100 -ofitval.dat

clean:
	-rm -f $(PROG) *.o core
tidy:
	-rm -f $(PROG) 

main.o : main.h ellipticity.h img.h besselfit.h psf.h normsinv.h ran2.h



