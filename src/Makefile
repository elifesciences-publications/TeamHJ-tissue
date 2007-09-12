all : ../bin/simulatorNew
#../bin/optimizer

# C++ Compiler
CXX = c++

BINDIR = ../bin

# Optimized
CXXFLAGS     = -O3 -DNDEBUG -funroll-loops -Wall -ansi -pedantic
LDFLAGS    = -O3 -DNDEBUG -funroll-loops -Wall -ansi -pedantic

#CXXFLAGS     = -O3 -funroll-loops -Wall -ansi -pedantic
#LDFLAGS    = -O3 -funroll-loops -Wall -ansi -pedantic

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
# CXXFLAGS     = -O0 -Wall -ansi -pedantic -g
# LDFLAGS    = -O0 -Wall -ansi -pedantic -g

# CXXFLAGS     = -O3 -Wall -ansi -pedantic -g
# LDFLAGS    = -O3 -Wall -ansi -pedantic -g

# List of object files.
OBJECTS = tissue.o cell.o wall.o vertex.o baseReaction.o growth.o mechanical.o mechanicalSpring.o baseCompartmentChange.o compartmentDivision.o compartmentRemoval.o network.o direction.o baseDirectionUpdate.o baseDirectionDivision.o directionUpdate.o directionDivision.o myConfig.o myFiles.o myRandom.o myTimes.o mySignal.o rungeKutta.o baseSolver.o math.o

# simulatorSolver

../bin/simulatorNew : $(OBJECTS) simulatorNew.o
	c++ simulatorNew.o $(OBJECTS) $(LDFLAGS) -o $(BINDIR)/simulatorNew

#../bin/optimizer: $(OBJECTS) optimizer.o
#	c++ optimizer.o $(OBJECTS) $(LDFLAGS) -o $(BINDIR)/optimizer


# Build rules

simulator.o : simulator.cc
#optimizer.o : optimizer.cc

$(OBJECTS): %.o: %.cc %.h

clean :
	rm *.o
