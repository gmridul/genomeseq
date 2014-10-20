CXX     = mpicxx
CXXFLAGS = -g -O3 -std=c++0x
LDLIBS = -lm -lrt
CLUSTALO = /home/mridul/Desktop/clustal-omega-1.2.1/src
CXXFLAGS += -I$(CLUSTALO)
#CFLAGS += -finstrument-functions
LDLIBS += -L$(CLUSTALO)/.libs -lclustalo -lgomp

OBJS =  cfg.o\
	cloneMates.o\
	dynamic.o\
	err.o\
	est.o\
	keys.o\
	master.o\
	mastercode.o\
	stree.o\
	suff.o\
	uFind.o\
	alignwrapper.o\

all: PaCE

pf: $(OBJS) pf.o
	gcc -c pf.c -o pf.o
	$(CC) $(OBJS) pf.o -o PaCE_pf $(LDLIBS)


PaCE: $(OBJS)
	$(CXX) $(OBJS) -o PaCE $(LDLIBS)

clean:
	rm -f *.o PaCE PaCE_pf
