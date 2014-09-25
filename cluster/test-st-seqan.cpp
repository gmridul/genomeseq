#include <iostream>
#include <seqan/index.h>
#include <vector>

using namespace seqan;

typedef Index< StringSet<CharString> , IndexEsa<> > TIndex;
typedef typename Iterator< TIndex, TopDown< ParentLinks<PostorderEmptyEdges> > >::Type TTreeIter;
typedef typename Infix< typename Fibre<TIndex, EsaBwt>::Type const >::Type	TOccsBWT;
typedef typename Infix< typename Fibre<TIndex, EsaSA>::Type const >::Type	TOccs;
typedef typename Iterator<TOccsBWT, Standard>::Type TIterBWT;
typedef typename Iterator<TOccs, Standard>::Type	TIter;
typedef typename Size<TIndex>::Type TSize;
typedef Pair<TSize> TPair;
typedef std::pair<int, int> MyPair;
typedef std::vector<MyPair> MyPairs;

template<class Sigma>
struct LSetInfo {
    Sigma  prev;
    int    skip;
    size_t read_id;
    size_t pos;
    LSetInfo() : prev('Z'), skip(0), read_id(0), pos(0) {}
};

typedef LSetInfo<char> TLSetInfo;

bool iter_compare(const TTreeIter & i, const TTreeIter & j) { 
    if (repLength(i) > repLength(j)) return true;
    if (repLength(i) == repLength(j)) return isLeaf(i);
    return false;
}

bool compare_lset(const TLSetInfo & i, const TLSetInfo & j) {
    if (i.prev < j.prev) return true;
    if ((i.prev == j.prev) && (i.read_id < j.read_id)) return true;
    if ((i.prev == j.prev) && (i.read_id == j.read_id) && (i.pos < j.pos)) return true;
    return false;
}

void print_lset(const std::vector<TLSetInfo>& lSet) {
    std::cout << "Printing lsets ...\n";
    std::cout << "  prev:";
    for (int i=0; i<lSet.size(); ++i) {
        std::cout << std::setw(3) << lSet[i].prev;
    }
    std::cout << "\n";
    std::cout << "  rdid:";
    for (int i=0; i<lSet.size(); ++i) {
        std::cout << std::setw(3) << lSet[i].read_id;
    }
    std::cout << "\n";
    std::cout << "  rpos:";
    for (int i=0; i<lSet.size(); ++i) {
        std::cout << std::setw(3) << lSet[i].pos;
    }
    std::cout << "\n";
}

void kway_merge(std::vector<TLSetInfo>& lSet, const std::vector<size_t>& positions) {
    std::cout << "kway merge: ";
    for (std::vector<size_t>::const_iterator i=positions.begin(); i!=positions.end(); ++i) {
        std::cout << *i << " ";
    }
    std::cout << "\n";
    std::sort(lSet.begin()+positions[0],lSet.begin()+positions[positions.size()-1], compare_lset); 
}

MyPairs & generate_pairs(std::vector<TLSetInfo>& lSet, int start, int stop) {
    //TODO signed (int) to unsiged (size_t)
    MyPairs pairs;
    int skip = 0;
    std::cout << "start=" << start << ", stop=" << stop << ", i=";
    for (int i=stop-2; i>=start; --i) {
        std::cout << " " << i << ", start=" << start;
        std::cout.flush();
        if (lSet[i].read_id == lSet[i+1].read_id) {
            lSet[i].skip = ++skip;
        } else {
            lSet[i].skip = skip = 0;
        }
        //print_lset(lSet);
    }
    std::cout << "\n";
    for (int i = start; i < stop-1; i += lSet[i].skip+1) {
        int j;
        for (j = lSet[i].skip+1; j < stop && (lSet[i].prev == lSet[j].prev); j += lSet[j].skip+1);
        for (; j < stop; j += lSet[j].skip+1) {
            pairs.push_back(std::make_pair(i,j));
        }
    }
}

int main ()
{
    StringSet<CharString> myStringSet;
    appendValue(myStringSet, "tobeornottobe");
//    appendValue(myStringSet, "tobeornottobe");
//    appendValue(myStringSet, "thebeeonthecomb");
//    appendValue(myStringSet, "beingjohnmalkovich");
    TIndex index(myStringSet);

    std::vector<TTreeIter> tree_iters;

    // Using Postorder iterator
    std::cout << "\nusing postorder";
    TTreeIter myIterator(index);

    // Top-down iterators start in the root node which is not the first node of a
    // postorder DFS. Thus we have to manually go to the DFS start with goBegin
    goBegin(myIterator);
    while (!atEnd(myIterator))
    {
		if (isLeaf(myIterator)) {
            std::cout << "leaf: ";
        } else {
            std::cout << "intl: ";
        }
        std::cout << ", " << representative(myIterator) << "(" << length(getOccurrences(myIterator)) << ")";
        std::cout << "[" << range(myIterator) << "]\n";
        tree_iters.push_back(myIterator);
        ++myIterator;
    }
    size_t numLeaves = length(index);
    std::cout << "num leaves = " << numLeaves << ", " << tree_iters.size() << "\n";
    std::vector<TLSetInfo> lSet;
    lSet.resize(numLeaves);

    //for (int i=0; i<tree_iters.size(); i++) {
    //    std::cout << representative(tree_iters[i]) << "\n";
    //}


    std::sort(tree_iters.begin(), tree_iters.end(), iter_compare); 

    indexRequire(index, EsaBwt());
//    typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

    for (int i=0; i<tree_iters.size(); i++) {
        TTreeIter& myIterator = tree_iters[i];
        TPair span = range(myIterator);
        TSize s = span.i1, t = span.i2;
        std::cout << representative(myIterator) << " " << range(myIterator) << " ";
        if (isLeaf(myIterator)) {
            std::cout << " leaf\n";
            TOccs occs = getOccurrences(myIterator);
            TOccsBWT bwts = getOccurrencesBwt(myIterator);
            TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());
            TIterBWT bw = begin(bwts, Standard());
            TIndex const &index = container(myIterator);
            typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

            while (oc != ocEnd) {
                lSet[s].pos = getSeqOffset(*oc);
                lSet[s].read_id = getSeqNo(*oc);
                if (posAtFirstLocal(*oc, limits)) {
                    lSet[s].prev = 'B';
                } else {
                    lSet[s].prev = *bw;
                }
                ++oc; ++bw; ++s;
            }
            std::sort(lSet.begin()+span.i1, lSet.begin()+span.i2, compare_lset);
        } else {
            std::cout << " int\n children are : \n  ";
            std::vector<size_t> positions;
            goDown(myIterator);
            do {
                std::cout << representative(myIterator) << " " << range(myIterator) << " ";
                positions.push_back(range(myIterator).i1);
            } while (goRight(myIterator));
            positions.push_back(range(myIterator).i2);
            std::cout << "\n";
            kway_merge(lSet, positions);
            generate_pairs(lSet, positions[0], positions[positions.size()-1]);
        }
        print_lset(lSet);
    }

    for (int i=0; i<numLeaves; i++) {
        std::cout << lSet[i].prev << " " << lSet[i].read_id << " " << lSet[i].pos << "\n";
    }



    return 0;
}

