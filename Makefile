CC = g++ 
COMPILER_FLAGS = -Wall -O3 -g

solver :
    $(CC) solver.cpp -o solver $(COMPILER_FLAGS)

.PHONY: solver clean

clean : 
    rm solver
