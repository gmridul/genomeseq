#ifndef __ALIGNWRAPPER_H__
#define __ALIGNWRAPPER_H__

#include "fastqwrapper.hpp"

extern "C" {
#include "clustal-omega.h"
}

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#define BOUND     2
#define THRESHOLD -30
#define KMERLEN   4
typedef seqan::String<char> TSequence;                     // sequence type
typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;   // align type

int match_seqs(const std::string &s1, const std::string &s2);

int myAlign(mseq_t *prMSeq, opts_t *prOpts);
void tree_to_align();
void testMuscleTree();
void testSetupSequences(mseq_t *prMSeq);
void setup_sequences(mseq_t *prMSeq, int nreads, Reads& reads);
void align_cluster(const Reads &reads, int num_leaves, uint *left, uint *right, uint *leaf_ids);

#endif // __ALIGNWRAPPER_H__ 

