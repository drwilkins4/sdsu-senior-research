COMPILER= gfortran

FLAGS = 

EXEC = program

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

math.o: math.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

physics.o: physics.f90 types.o math.o
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 types.o physics.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 types.o read_write.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)