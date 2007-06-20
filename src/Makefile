all : ../bin/simulator
#../bin/optimizer

# C++ Compiler
CXX = c++

BINDIR = ../bin

# Optimized
CXXFLAGS     = -O3 -DNDEBUG -funroll-loops -Wall -ansi -pedantic
LDFLAGS    = -O3 -DNDEBUG -funroll-loops -Wall -ansi -pedantic

# Including all assert calls
#CXXFLAGS     = -O3
#LDFLAGS    = -O3

# For compiling code for the cluster
#CXXFLAGS     = -O3
#LDFLAGS    = -O3 -static -all-static

# for gprof time testing
#CXXFLAGS     = -O3 -pg -DNDEBUG
#LDFLAGS    = -O3 -pg -DNDEBUG

#CPPPATH      = ~/include
#LIBPATH      = ~/lib

# For development 
#CXXFLAGS     = -Wall -ansi -pedantic -g
#LDFLAGS    = -Wall -ansi -pedantic -g

# List of object files.
OBJECTS = tissue.o cell.o wall.o vertex.o baseReaction.o growth.o mechanical.o baseCompartmentChange.o compartmentDivision.o compartmentRemoval.o network.o direction.o baseDirectionUpdate.o baseDirectionDivision.o directionUpdate.o directionDivision.o myConfig.o myFiles.o myRandom.o

# simulatorSolver

../bin/simulator : $(OBJECTS) simulator.o
	c++ simulator.o $(OBJECTS) $(LDFLAGS) -o $(BINDIR)/simulator

#../bin/optimizer: $(OBJECTS) optimizer.o
#	c++ optimizer.o $(OBJECTS) $(LDFLAGS) -o $(BINDIR)/optimizer


# Build rules

simulator.o : simulator.cc
#optimizer.o : optimizer.cc

$(OBJECTS): %.o: %.cc %.h

clean :
	rm *.o
