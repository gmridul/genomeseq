#include <iostream>
#include <seqan/index.h>
#include <vector>
#include <algorithm>
#include <utility>

using namespace seqan;
using namespace std;

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
    size_t read_id;
    size_t pos;
    LSetInfo() : prev('Z'), read_id(0), pos(0) {}
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

void print_lset(const std::vector<TLSetInfo>& lSet, const string msg="lset") {
    std::cout << "Printing " << msg << " ...\n";
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


struct qelement {
    int pos;
    int part;
    qelement(int _pos, int _part) : pos(_pos), part(_part) {}
};

class mycomparison {
    public:
        mycomparison(vector<TLSetInfo>& lSet): _lSet(lSet) {}
        bool operator() (const qelement &lhs, const qelement& rhs) {
            return ((_lSet[lhs.pos].prev>_lSet[rhs.pos].prev) || ((_lSet[lhs.pos].prev==_lSet[rhs.pos].prev) && (_lSet[lhs.pos].read_id>_lSet[rhs.pos].read_id)));
        }
    private:
        vector<TLSetInfo> & _lSet;
};

void kway_merge(std::vector<TLSetInfo>& lSet, std::vector<TLSetInfo>& lSetCopy, std::vector<size_t>& pos) {
    std::cout << "kway merge: ";
    for (std::vector<size_t>::const_iterator i=pos.begin(); i!=pos.end(); ++i) {
        std::cout << *i << " ";
    }
    std::cout << "\n";
    //    std::sort(lSet.begin()+positions[0],lSet.begin()+positions[positions.size()-1], compare_lset); 

    int n = pos.size();
    vector<TLSetInfo> temp(lSet.begin()+pos[0],lSet.begin()+pos[n-1]);
    print_lset(temp, "temp");


    priority_queue<qelement,vector<qelement>,mycomparison> outQueue(temp);

    for (int i=0; i<n-1; i++) {
        outQueue.push(qelement(pos[i]-pos[0], i));
    }

    int s=pos[0];
    while (!outQueue.empty()) {
        const qelement& lowest = outQueue.top();
        //std::cout << "lowest.pos=" << outQueue.top().pos << ", lowest.part=" << outQueue.top().part << "\n";
        lSet[s++]=temp[lowest.pos];
        if(lowest.pos+pos[0]+1 < pos[lowest.part+1]) {
            qelement qe(lowest.pos+1, lowest.part);
            //td::cout << "subsequent pushing (" << qe.pos << "," << qe.part << ")\n";
            outQueue.pop();
            outQueue.push(qe);
        } else {
            outQueue.pop();
        }
        //print_lset(temp, "temp");
        //print_lset(lSet);
    }
}

void generate_pairs_from_two(MyPairs& pairs, const std::vector<TLSetInfo>& lSet, size_t start1, size_t end1, size_t start2, size_t end2) {
    for (size_t i = start1; i < end1; ++i) {
        for (size_t j = start2; j < end2; ++j) {
            std::cout << "i=" << i << ", j=" << j << "\n";
            if ((lSet[i].read_id != lSet[j].read_id) &&
                ((lSet[i].prev < lSet[j].prev) || ((lSet[i].prev == 'B') && (lSet[j].prev == 'B')))) {
                pairs.push_back(std::make_pair(i,j));
            }
        }
    }
}

MyPairs & generate_pairs(const std::vector<TLSetInfo>& lSet, const std::vector<size_t>& positions) {
    MyPairs pairs;
    for (size_t i=0; i<positions.size()-2; ++i) {
        size_t start1 = positions[i], end1 = positions[i+1];
        for (size_t j=i+1; j<positions.size()-1; ++j) {
            size_t start2 = positions[j], end2 = positions[j+1];
            generate_pairs_from_two(pairs, lSet, start1, end1, start2, end2);
        }
    }

    std::cout << "all pairs:";
    for (int i=0; i<pairs.size(); ++i) {
        int f = pairs[i].first, s = pairs[i].second;
        std::cout << "<" << lSet[f].read_id << ":" << lSet[f].pos << "," << lSet[s].read_id << ":" << lSet[s].pos << "> ";
    }
    std::cout << "\n";
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
    std::vector<TLSetInfo> lSet, lSetCopy;
    lSet.resize(numLeaves);
    lSetCopy.resize(numLeaves);

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
        std::vector<size_t> positions;
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
            goDown(myIterator);
            do {
                std::cout << representative(myIterator) << " " << range(myIterator) << " ";
                positions.push_back(range(myIterator).i1);
            } while (goRight(myIterator));
            positions.push_back(range(myIterator).i2);
            std::cout << "\n";
            //generate_pairs(lSet, positions);
            kway_merge(lSet, lSetCopy, positions);
        }
        print_lset(lSet);
    }

    for (int i=0; i<numLeaves; i++) {
        std::cout << lSet[i].prev << " " << lSet[i].read_id << " " << lSet[i].pos << "\n";
    }



    return 0;
}

