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
DOS.e: Diag.o main.o
	$(CC) $(FLAGS) Diag.o main.o $(LIBS) -o DOS.e

Diag.o: Diag.f90 
	$(CC) $(FLAGS) -c Diag.f90

main.o: main.f90 Diag.o
	$(CC) $(FLAGS) -c main.f90

clean:
	rm -f *.o *.mid *~

endall:
	rm *.text
