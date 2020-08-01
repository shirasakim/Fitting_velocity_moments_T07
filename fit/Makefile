## Compiler, tools and options
CC 	= gcc
CFLAGS 	= -O3 -Wall #-I/home/shirskms/usr_gnu/local/include
INCDIRS = #-I/opt/local/include
LIBDIRS = #-L/home/shirskms/usr_gnu/local/lib #-L/opt/local/lib
LIBS 	= -lm -lgsl -lgslcblas

SRCH = #my_cpgplot.h
fit_mean= fit_mean.o
fit_variance= fit_variance.o

fit_mean: $(fit_mean)
	$(CC) $(fit_mean) -o fit_mean $(CFLAGS) $(LIBDIRS) $(LIBS)

fit_variance: $(fit_variance)
	$(CC) $(fit_variance) -o fit_variance $(CFLAGS) $(LIBDIRS) $(LIBS)

## Implicit rules

.SUFFIXES: .c .o

## Build rules ##

## Compile

.c.o: $*.c
	$(CC) -c $*.c $(CFLAGS) $(INCDIRS)

## Clean

clean:
	rm -f $(ALL) $(TARGET) *.o *~ fit_mean fit_variance
