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
#define BOUND 2
#define THRESHOLD -30
typedef seqan::String<char> TSequence;                     // sequence type
typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;   // align type

int myAlign(mseq_t *prMSeq, opts_t *prOpts);
void tree_to_align();
void testMuscleTree();
void testSetupSequences(mseq_t *prMSeq);
void setup_sequences(mseq_t *prMSeq, int nreads, Reads& reads);
void align_cluster(int n, int ni, uint *left, uint *right, float *leftLength, float *rightLength, uint *leafIds, char **leafNames, uint root, mseq_t *prMSeq);

#endif // __ALIGNWRAPPER_H__ 

