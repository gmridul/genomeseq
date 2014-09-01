#ifndef __FASTQWRAPPER_H__
#define __FASTQWRAPPER_H__

#include <unordered_map>
#include <vector>
#include <string>

class llist {
    public:
        llist* next;
        char side;
        int pos;
        int entrynum;
        llist* cluster;
        llist() {
            side='c';
        }
};
typedef std::unordered_map<int64_t,llist* > HashTable;
typedef std::vector<std::string> Reads;

int read_fastq(char *fname1, char *fname2, Reads &left_reads, Reads &right_reads, int k, HashTable &hashtab);
void print_hashtab(const HashTable &hashtab);

#endif // __FASTQWRAPPER_H__ 

