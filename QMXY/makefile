#
CC=gcc -O2 -Wall -I/usr/local/include/gsl
CD=gcc -o 
CFLAGS= -L/usr/local/lib -lm -lgsl -lgslcblas
EXTRAS=
FFLAGS= -c
# define object files

OBJ1 = qmxyv0.o random.o 
#-------------------------------------------------------------------

$(OBJ1) : $(HFILES)

qmxy: $(OBJ1)
	$(CD) qmxy $(OBJ1) $(EXTRAS) $(CFLAGS)

# --------------------------------------------------------------------
clean:
	rm $(OBJ1) 

# declare dependencies
.c.o: $(HFILES)
	$(CC) $(FFLAGS) $*.c
