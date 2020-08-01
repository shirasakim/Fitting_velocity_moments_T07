## Compiler, tools and options
CC 	= gcc
CFLAGS 	= -O3 -Wall #-I/home/shirskms/usr_gnu/local/include
INCDIRS = #-I/opt/local/include
LIBDIRS = #-L/home/shirskms/usr_gnu/local/lib #-L/opt/local/lib
LIBS 	= -lm

SRCH = #my_cpgplot.h
gen_pdf= gen_pdf_rockstar.o

gen_pdf_rockstar: $(gen_pdf)
	$(CC) $(gen_pdf) -o gen_pdf_rockstar $(CFLAGS) $(LIBDIRS) $(LIBS)

## Implicit rules

.SUFFIXES: .c .o

## Build rules ##

## Compile

.c.o: $*.c
	$(CC) -c $*.c $(CFLAGS) $(INCDIRS)

## Clean

clean:
	rm -f $(ALL) $(TARGET) *.o *~ gen_pdf_rockstar
