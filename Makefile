CXX = g++
INC_DIR = ../nlohmann
CXXFLAGS = -std=c++11 -O1 -Wall -pedantic -g -fopenmp -I "/usr/include/eigen3" -I$(INC_DIR)
# Settings for UON HPC
#CXXFLAGS = -std=gnu++11 -O2 -Wall -Wextra -pedantic -I /cm/local/apps/boost/1.65.1/include
#CXXFLAGS = -std=c++11 -O1 -Wall -pedantic -g -fopenmp -I "/cm/software/apps/eigen/gcc8/3.3.7/include/Eigen/" -I$(INC_DIR)
LDFLAGS =  -lgomp #-pthread

LIST = main data eval pop agent cont_frac pop_init \
       pop_evolve pop_recomb pop_mutate local_search linear objective \
       debug helper args 

SRC = $(addsuffix .cpp, $(LIST))
OBJ = $(addprefix bin/, $(addsuffix .o, $(LIST)))

all: main

bin/%.o : %.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

main: $(OBJ)
	$(CXX) $^ $(LDFLAGS) -o bin/main

clean:
	rm -f bin/*

.PHONY : all clean
