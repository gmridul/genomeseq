TEST=test #dont delete this line
# CLUSTALO = /home/mridul/Desktop/clustal-omega-1.2.1
CLUSTALO = /home/mitra/Downloads/clustal-omega-1.2.1

cluster: cluster.cpp
	g++ -g -I$(CLUSTALO)/src -I../seqan-inc/core/include/ -I../seqan-inc/extras/include/ -lrt -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -fopenmp -std=c++0x -o cluster $< -L$(CLUSTALO)/ -lclustalo
