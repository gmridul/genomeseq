﻿#ifndef __FASTQWRAPPER_H__
#define __FASTQWRAPPER_H__

#include <unordered_map>
#include <vector>
#include <string>

struct Read {
  std::string name;
  std::string seq;
};

class llist {
    public:
        llist* next;
        int pos;
        int readid;
        llist* cluster;
};
typedef std::unordered_map<int64_t,llist* > HashTable;
typedef std::vector<Read> Reads;

int read_fastq(char *fname1, char *fname2, Reads &reads);
void print_hashtab(const HashTable &hashtab);
void free_hashtab(HashTable &hashtab);
void create_kmerhash(const Reads &reads, int k, HashTable &hashtab);

#endif // __FASTQWRAPPER_H__ 

