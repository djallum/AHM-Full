# define the compiler to use 
CC = gfortran

# compiler flags:
#      -g -fcheck=all -Wall is for debugging
#      -O0 and -O3 defines the level of optimization (0 is low 3 is high) not sure what this is for
FLAGS =  -g -fcheck=all -Wall
#FLAGS = -O3

# define libraries
LIBS = -llapack

# Define date
DD = 180330
DOS.e: $(DD)Diag.o $(DD)main.o
	$(CC) $(FLAGS) $(DD)Diag.o $(DD)main.o $(LIBS) -o DOS.e

$(DD)Diag.o: $(DD)Diag.f90 
	$(CC) $(FLAGS) -c $(DD)Diag.f90

$(DD)main.o: $(DD)main.f90 $(DD)Diag.o
	$(CC) $(FLAGS) -c $(DD)main.f90

clean:
	rm -f *.o *.mid *~
