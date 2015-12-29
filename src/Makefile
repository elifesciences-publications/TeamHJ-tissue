# Makefile copied from organism.

# This Makefile should automatically keeps track of dependencies in the
# code. It saves dependencies in *.d files and includes them whenever
# necessary. The dependency part is based on Scott McPeak's Makefile at
# http://www.cs.berkeley.edu/~smcpeak/autodepend/autodepend.html. The
# other parts are based on Simon's Makefile and tons of obscure
# webpages I no longer recall.  
#
#'make' build executable files
# 	'bin/simulator' and 'bin/optimizer' 
#
#'make debug'	Compiles with -g and no optimization 
#
#'make gprof'	Compiles with -pg and no optimization 
#
#'make clean' 	removes all .o, executables and .d 

SHELL = /bin/sh
#Compilers and flags
CXX = g++
CXXFLAGS = -g -O3 -DNDEBUG -Wall -pedantic -ansi -std=c++11
LDFLAGS = -g -O3 -DNDEBUG -Wall -pedantic -ansi -std=c++11
#CXXFLAGS = -g -O3 -DNDEBUG -Wall -pedantic   -std=c++11 -stdlib=libc++ -I/sw/include/
#LDFLAGS = -g -O3 -DNDEBUG -Wall -pedantic 
#CXXFLAGS = -g -Wall -pedantic -ansi
#LDFLAGS = -g -Wall -pedantic -ansi 
#
# uncomment two rows below for MAC OS X
#
#CXXFLAGS = -g -O3 -DNDEBUG -Wall -pedantic -std=c++11 -I/sw/include/ -I/opt/local/include/
#LDFLAGS = -g -O3 -DNDEBUG -Wall -pedantic

#Sourcefiles etc.
SIM_SRC = simulator/simulator.cc
SIM_OBJ = $(SIM_SRC:.cc=.o)
CONV_SRC = tools/converter.cc
CONV_OBJ = $(CONV_SRC:.cc=.o)
SRC_DIR = ./ ./ply/
SRCS := $(wildcard $(SRC_DIR:/=/*.cc))
OBJS = $(SRCS:.cc=.o)
DEP = makedep.rule

#Binaries
SIMULATOR = ../bin/simulator
CONVERTER = ../bin/converter

all: $(SIMULATOR) $(CONVERTER)

$(SIMULATOR): $(OBJS) $(SIM_OBJ)
	$(CXX) $(SIM_OBJ) $(OBJS) $(LDFLAGS) -o $(SIMULATOR) 

$(CONVERTER): $(OBJS) $(CONV_OBJ)
	$(CXX) $(CONV_OBJ) $(OBJS) $(LDFLAGS) -o $(CONVERTER) 

$(SIM_OBJ) : $(SIM_SRC)

$(CONV_OBJ) : $(CONV_SRC)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# Compile and generate dependency info.  The two first lines fixes a
# bug in gcc(?), and prints dir/foo.o: dir/foo.cc... in the dependency
# file instead of foo.o: dir/foo.cc... as is done by gcc. The other
# two lines use a trick to get make to detect name changes etc.
%.o: %.cc
	$(CXX) -o $*.o -c $(CXXFLAGS) $*.cc
	$(CXX) -MM $(CXXFLAGS) $*.cc > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

debug:
	make "CXXFLAGS = -g -pedantic -Wall -ansi" "LDFLAGS = -g -pedantic -Wall -ansi" 		
gprof:
	make "CXXFLAGS = -pg -pedantic -Wall -ansi" "LDFLAGS = -pg -pedantic -Wall -ansi" 		


.PHONY: clean

clean:
	rm -f $(OBJS)
	rm -f $(SIM_OBJ)
	rm -f $(CONV_OBJ)
	rm -f $(SIMULATOR)
	rm -f $(CONVERTER)
	rm -f *.d
	rm -f */*.d
