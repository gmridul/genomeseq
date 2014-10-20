#ifndef __ALIGNWRAPPER_H__
#define __ALIGNWRAPPER_H__

#include <vector>
#include <string>
extern "C" {
#include "clustal-omega.h"
}
struct Read {
  std::string name;
  std::string seq;
};

typedef std::vector<Read> Reads;

void setup_sequences(mseq_t *prMSeq, int nreads, Reads& reads);
void align_cluster(const Reads &reads, int num_leaves, uint *left, uint *right, uint *leaf_ids);

#endif // __ALIGNWRAPPER_H__ 

