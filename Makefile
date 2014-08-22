cluster: cluster.cpp
	g++ -g -I/home/mridul/Desktop/clustal-omega-1.2.1/src -I./seqan-inc/core/include/ -I./seqan-inc/extras/include/ -lrt -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -fopenmp -std=c++0x -o cluster $< -L/usr/local/lib -lclustalo
