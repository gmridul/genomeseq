cluster: cluster.cpp
	g++ -g -I./seqan-inc/core/include/ -I./seqan-inc/extras/include/ -lrt -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -std=c++0x -o cluster $<
