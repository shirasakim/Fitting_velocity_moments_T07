
GSLINC = /usr/local/include
GSLLIB = /usr/local/lib

# directory where fftw3.h exists
FFTINC = /usr/local/include
# directory where libcfitsio.a exists
FFTLIB = /usr/local/lib

OBJ1 = save_mass_variance_EH.o\
	nrutil.o qromb.o polint.o trapzd.o spline.o
EXE1 = save_mass_variance_EH

#linker
LD = g++ 

#c compiler
CC = g++

#cpp compiler
CPP = g++

#CPP compiler options
CPPOPT = -O3 #

#C compiler options
COPT = -O3 -s # -bmaxstack:512000000 -g

#linking
LIBS    = -lgsl -lgslcblas -lm
#compilation
.C.o:
	$(CC) $(COPT) $(FCFLAGS) $(CPPFLAGS) -c $.C

%.o: %.C
	$(CPP) $(CPPOPT) -c $*.C

$(EXE1): $(OBJ1)
	$(LD) $(CPPOPT) -o $(EXE1) $(OBJ1) $(LIBS)


.FAILED:
.DONE:

clean:
	rm -f *.o *~ $(EXE1)

#dependencies

