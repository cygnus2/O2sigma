#
CC=gcc -O2 -Wall -I/home/debasish/dev/include
CD=gcc -o 
CFLAGS= -L/home/debasish/dev/lib -lm -lgsl -lgslcblas
EXTRAS=
FFLAGS= -c
# define object files

OBJ1 = xy4dobc.o random.o 
#-------------------------------------------------------------------

$(OBJ1) : $(HFILES)

xy4dobc: $(OBJ1)
	$(CD) xy4dobc $(OBJ1) $(EXTRAS) $(CFLAGS)

# --------------------------------------------------------------------
clean:
	rm $(OBJ1) 

# declare dependencies
.c.o: $(HFILES)
	$(CC) $(FFLAGS) $*.c
