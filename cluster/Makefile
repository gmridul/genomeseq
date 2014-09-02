CXX     = g++
CXXFLAGS = -g -Wall -Wno-long-long -pedantic -Wno-variadic-macros -fopenmp -std=c++0x -DNDEBUG # -O3
LDLIBS = -lm -lrt 

# CLUSTALO = /home/mridul/Desktop/clustal-omega-1.2.1
CLUSTALO = /home/mitra/Downloads/clustal-omega-1.2.1/src
CXXFLAGS += -I$(CLUSTALO) -I../seqan-inc/core/include/ -I../seqan-inc/extras/include/
LDLIBS += -L$(CLUSTALO)/.libs -lclustalo -lgomp

OBJS =  cluster.o\
	fastqwrapper.o\
	alignwrapper.o

all: cluster

cluster: $(OBJS)
	$(CXX) $(OBJS) -o cluster $(LDLIBS)

clean:
	rm -f *.o cluster

