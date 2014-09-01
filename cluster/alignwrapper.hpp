#ifndef __ALIGNWRAPPER_H__
#define __ALIGNWRAPPER_H__


extern "C" {
#include "clustal-omega.h"
}

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#define BOUND 2
#define THRESHOLD -10
typedef seqan::String<char> TSequence;                     // sequence type
typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;   // align type

int myAlign(mseq_t *prMSeq, opts_t *prOpts);
void tree_to_align(char* t);

#endif // __ALIGNWRAPPER_H__ 

