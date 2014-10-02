
#include <iostream>

#include "fastqwrapper.hpp"
#include "alignwrapper.hpp"
#include "dswrapper.hpp"

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
    //std::cout << "kway merge: ";
    //for (std::vector<size_t>::const_iterator i=pos.begin(); i!=pos.end(); ++i) {
    //    std::cout << *i << " ";
    //}
    //std::cout << "\n";
    //    std::sort(lSet.begin()+positions[0],lSet.begin()+positions[positions.size()-1], compare_lset); 

    int n = pos.size();
    vector<TLSetInfo> temp(lSet.begin()+pos[0],lSet.begin()+pos[n-1]);
    //print_lset(temp, "temp");


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


struct BinaryTreeInfo {
    uint root;
    uint num_leaves;
    uint *left;
    uint *right;
    uint *leafIds;
    void print() {
        std::cout << "Number of leaves : " << num_leaves << "\n";
        std::cout << "LeafIds: ";
        for (uint j=0; j<num_leaves; ++j) {
            std::cout << std::setw(5) << leafIds[j] << " ";
        }
        std::cout << "\n";
        std::cout << "Left   : ";
        for (uint j=0; j<num_leaves-1; ++j) {
            std::cout << std::setw(5) << left[j] << " ";
        }
        std::cout << "\n";
        std::cout << "Right  : ";
        for (uint j=0; j<num_leaves-1; ++j) {
            std::cout << std::setw(5) << right[j] << " ";
        }
        std::cout << "\n";
    }
};

class BinaryForest {
    public:
        BinaryForest(int num_leaves): _internal_count_forest(0) {
            _forest.num_leaves = num_leaves;
            _forest.left = new uint[num_leaves-1];
            _forest.right = new uint[num_leaves-1];
            _forest.leafIds = new uint[num_leaves];
            for (auto i = 0; i < num_leaves-1; ++i) {
                _forest.left[i] = _forest.right[i] = num_leaves+i;
            }
            for (auto i = 0; i < num_leaves; ++i) {
                _forest.leafIds[i] = i;
            }
            _cur_tree.num_leaves = 0;
            _cur_tree.left = new uint[num_leaves-1];
            _cur_tree.right = new uint[num_leaves-1];
            _cur_tree.leafIds = new uint[num_leaves];
        }
        ~BinaryForest() {
            delete [] _forest.left;
            delete [] _forest.right;
            delete [] _cur_tree.left;
            delete [] _cur_tree.right;
            delete [] _cur_tree.leafIds;
        }
        uint addNode(uint lnode, uint rnode) {
            _forest.left[_internal_count_forest] = lnode;
            _forest.right[_internal_count_forest] = rnode;
            return _forest.num_leaves + _internal_count_forest++;
        }
        BinaryTreeInfo* getTree(uint root) {
            _internal_count_tree = 0;
            _cur_tree.num_leaves = 0;
            extractTree(root);
            adjustArray(_cur_tree.left, _forest.num_leaves, _cur_tree.num_leaves);
            adjustArray(_cur_tree.right, _forest.num_leaves, _cur_tree.num_leaves);
            return &_cur_tree;
        }

    private:
        void adjustArray(uint *array, uint old, uint curr) {
            for (uint i=0; i<curr-1; ++i) {
                if (array[i] >= old) array[i] += curr - old;
            }
        }

        uint extractTree(uint node) {
            //      std::cout << "extract tree with node " << node << "\n";
            uint current_id;
            if (node < _forest.num_leaves) {
                current_id = _cur_tree.num_leaves++;
                _cur_tree.leafIds[current_id] = node;
                return current_id;
            }
            current_id = _internal_count_tree++;
            node -= _forest.num_leaves;
            _cur_tree.left[current_id] = extractTree(_forest.left[node]);
            _cur_tree.right[current_id] = extractTree(_forest.right[node]);
            return _forest.num_leaves + current_id;
        }

        uint _internal_count_forest;
        uint _internal_count_tree;
        BinaryTreeInfo _forest;
        BinaryTreeInfo _cur_tree;
};

std::string reverse_complement(std::string read) {
    int len = read.length();
    std::string reverse;
    reverse.resize(len);
    int rlen=0;
    while(len-->0) {
        if(read[len]=='A') reverse[rlen]='T';
        else if(read[len]=='C') reverse[rlen]='G';
        else if(read[len]=='T') reverse[rlen]='A';
        else if(read[len]=='G') reverse[rlen]='C';
        rlen++;
    }
    return reverse;
}


#define RANK(x) (x->readid)

class Clusters {
    public:
        Clusters(int n) : num_elements(n), forest(num_elements), myrank(elements), parent(elements), ds(&myrank, &parent) {
            set_to_node = new uint[num_elements];
            elements.reserve(num_elements);
            for (size_t i = 0; i < elements.capacity(); ++i) {
                elements.push_back(Element(i));
            }
            for (size_t i = 0; i < elements.size(); ++i) {
                elements[i].dsID = i;
                set_to_node[i] = i;
            }
            for (size_t i = 0; i < elements.size(); ++i) {
                ds.make_set(elements.at(i));
            }
        }
        ~Clusters() {
            delete [] set_to_node;
        }
        void unite_by_parent(Element& xp, Element& yp) {
            int nodex = set_to_node[xp.dsID];
            int nodey = set_to_node[yp.dsID];
            ds.link(xp, yp);
            Element zp = ds.find_set(elements[xp.dsID]);
            set_to_node[zp.dsID] = forest.addNode(nodex, nodey);
        }
        void unite(int x, int y) {
            Element xp = ds.find_set(elements[x]);
            Element yp = ds.find_set(elements[y]);
            if (xp != yp) {
                unite_by_parent(xp, yp);
            }
        }
        void finalize_clusters() {
            ds.compress_sets(elements.begin(), elements.end());
            std::sort(elements.begin(), elements.end(), compareByParent);
            size_t first = 0;
            while (first < elements.size()) {
                size_t currentParent = elements.at(first).dsParent;
                size_t last = first;
                while (last < elements.size() && elements.at(last).dsParent == currentParent) {
                    ++last;
                }
                cluster_sizes.push_back(last-first);
                first = last;
            }
        }
        void print_elements() {
            printElements(elements);
        }
        Element & getElement(int id) {
            return elements[id];
        }
        /*void cluster_reads(const Reads &reads, HashTable &hashtab) {
            std::cout << "doing clusteing \n";
            int **checked = new int*[num_elements];
            for (int i=0; i<num_elements; i++) {
                checked[i] = new int[num_elements];
                for (int j=0; j<num_elements; j++) {
                    checked[i][j]=0;
                }
            }
            for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
                //std::cout << "checking hash table  " << it->first << "  " << it->second->readid << "\n";
                for (auto x = it->second; (x != NULL) && (x->next != NULL); x = x->next) {
                    int xrank = RANK(x);
                    for (auto y = x->next; y != NULL; y = y->next) {
                        int yrank = RANK(y);
                        //std::cout << "x=" << xrank << "(" << (xrank & ~0x1) << "), y=" << yrank << "(" << (yrank & ~0x1) << ")\n";
                        //                    if ((xrank & ~0x1) == (yrank & ~0x1)) {
                        //                        // from same pair
                        //                        continue;
                        //                    }
                        if(checked[xrank][yrank]) continue;
                        checked[xrank][yrank] = checked[yrank][xrank] = 1;

                        Element xp = ds.find_set(elements[xrank]);
                        Element yp = ds.find_set(elements[yrank]);
                        //Element xp_pair = ds.find_set(elements[xrank ^ 0x1]);
                        //Element yp_pair = ds.find_set(elements[yrank ^ 0x1]);


                        //if (xp.dsID == yp.dsID || xp.dsID == yp_pair.dsID || xp_pair.dsID == yp.dsID) continue;
                        if (xp.dsID == yp.dsID) continue;

                        //std::cout << reads[yrank] << "\n" << reverse_complement(reads[yrank]) << "\n";
                        if (match_seqs(reads[xrank].seq, reads[yrank].seq) >= THRESHOLD) {
                            unite_by_parent(xp, yp);
                        }
                        //                    else if (match_seqs(reads[xrank], reverse_complement(reads[yrank])) >= THRESHOLD) {
                        //                        unite_by_parent(xp, yp);
                        //                    }
                    }
                }
            }
            std::cout << ".... done \n";
            for (int i=0; i<num_elements; i++)
                delete [] checked[i];
            delete [] checked;
        }*/

        void generate_pairs_from_two(const Reads &reads, const std::vector<TLSetInfo>& lSet, size_t start1, size_t end1, size_t start2, size_t end2) {
            for (size_t i = start1; i < end1; ++i) {
                for (size_t j = start2; j < end2; ++j) {
                 //   std::cout << "i=" << i << ", j=" << j << "\n";
                    if ((lSet[i].read_id != lSet[j].read_id) && ((lSet[i].prev < lSet[j].prev) || ((lSet[i].prev == 'B') && (lSet[j].prev == 'B')))) {
                        Element xp = ds.find_set(elements[lSet[i].read_id]);
                        Element yp = ds.find_set(elements[lSet[j].read_id]);
                        if (xp.dsID == yp.dsID) return;
                        if (match_seqs(reads[lSet[i].read_id].seq, reads[lSet[j].read_id].seq) >= THRESHOLD) {
                            unite_by_parent(xp, yp);
                        }
                    }
                }
            }
        }

        void generate_pairs(const Reads &reads, const std::vector<TLSetInfo>& lSet, const std::vector<size_t>& positions) {
            for (size_t i=0; i<positions.size()-2; ++i) {
                size_t start1 = positions[i], end1 = positions[i+1];
                for (size_t j=i+1; j<positions.size()-1; ++j) {
                    size_t start2 = positions[j], end2 = positions[j+1];
                    generate_pairs_from_two(reads, lSet, start1, end1, start2, end2);
                }
            }
        }


        void cluster_reads_st(const Reads &reads) {
            std::cout << "doing clusteing \n";

            StringSet<CharString> myStringSet;
            for (int i=0; i<reads.size(); i++) {
                appendValue(myStringSet, reads[i].seq);
            }
            TIndex index(myStringSet);

            std::vector<TTreeIter> tree_iters;

            // Using Postorder iterator
            std::cout << "\nusing postorder";
            TTreeIter myIterator(index);
            goBegin(myIterator);
            while (!atEnd(myIterator)) {
                if (repLength(myIterator) >= KMERLEN) {
                    tree_iters.push_back(myIterator);
                }
                ++myIterator;
            }
            size_t numLeaves = length(index);
            //std::cout << "num leaves = " << numLeaves << ", " << tree_iters.size() << "\n";
            std::vector<TLSetInfo> lSet, lSetCopy;
            lSet.resize(numLeaves);
            lSetCopy.resize(numLeaves);

            // TODO: use bucket sort or something
            std::sort(tree_iters.begin(), tree_iters.end(), iter_compare); 

            indexRequire(index, EsaBwt());

            for (int i=0; i<tree_iters.size(); i++) {
                TTreeIter& myIterator = tree_iters[i];
                TPair span = range(myIterator);
                TSize s = span.i1, t = span.i2;
                //std::cout << representative(myIterator) << " " << range(myIterator) << " ";
                std::vector<size_t> positions;
                if (isLeaf(myIterator)) {
                    //std::cout << " leaf\n";
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
                    //std::cout << " int\n children are : \n  ";
                    goDown(myIterator);
                    do {
                        //std::cout << representative(myIterator) << " " << range(myIterator) << " ";
                        positions.push_back(range(myIterator).i1);
                    } while (goRight(myIterator));
                    positions.push_back(range(myIterator).i2);
                    //std::cout << "\n";
                    generate_pairs(reads, lSet, positions);
                    kway_merge(lSet, lSetCopy, positions);
                }
                print_lset(lSet);
            }

            //for (int i=0; i<numLeaves; i++) {
            //    std::cout << lSet[i].prev << " " << lSet[i].read_id << " " << lSet[i].pos << "\n";
            //}
        }


        std::vector<int>& get_cluster_sizes() { return cluster_sizes; }
        BinaryTreeInfo *getTree(int id) { return forest.getTree(set_to_node[elements[id].dsParent]); }

    private:
        int num_elements;
        BinaryForest forest;
        MyRank myrank;
        Parent parent;
        DisjointSets ds;
        std::vector<Element> elements;
        std::vector<int> cluster_sizes;
        uint *set_to_node;
};

//#define TEST_SIMPLE

int main(int argc, char*  argv[]) {

    Reads reads;

#ifdef TEST_SIMPLE
    reads.push_back("AAA");
    reads.push_back("AC");
    reads.push_back("AAT");
    reads.push_back("AAG");
    reads.push_back("AAAT");
    reads.push_back("ACAT");
    reads.push_back("AAT");
    reads.push_back("ATAA");
    reads.push_back("AATA");
#else
    read_fastq(argv[1], argv[2], reads);
    HashTable hashtab;
    create_kmerhash(reads, KMERLEN, hashtab);
    //print_hashtab(hashtab);
#endif

    //for (int i=0; i<reads.size(); i++) {
    //    std::cout << reads[i] << "\n";
    //}
    int n = reads.size();
    Clusters cls(n);

#ifdef TEST_SIMPLE
    cls.unite(0, 1);
    cls.unite(2, 3);
    cls.unite(2, 4);
    cls.unite(1, 4);
    cls.unite(5, 6);
    cls.unite(7, 8);
    cls.unite(5, 8);
#else
    cls.cluster_reads_st(reads);
    //free_hashtab(hashtab);
#endif

    cls.finalize_clusters();
    //cls.print_elements();
    std::string fcls_name("cluster_info.dat");
    std::ofstream fcls(fcls_name);
    if (!fcls.is_open()) {
        std::cerr << "ERROR: Could not open file " << fcls_name << " for writing.\n"; 
    }
    std::vector<int> sizes = cls.get_cluster_sizes();
    int id=0;

    for (size_t i = 0; i < sizes.size(); ++i) {
        std::cout << "Cluster " << i << " has " << sizes[i] << " elements:";
        for (int j=0; j<sizes[i]; ++j) {
            std::cout << " " << cls.getElement(id+j).dsID;
        }
        std::cout << "\n";
        fcls << reads[cls.getElement(id).dsID].name;
        for (int j=1; j<sizes[i]; ++j) {
            fcls << " " << reads[cls.getElement(id+j).dsID].name;
        }
        fcls << "\n";
        if (sizes[i] > 1) {
             BinaryTreeInfo *tree = cls.getTree(id);
             //tree->print();
             align_cluster(reads, sizes[i], tree->left, tree->right, tree->leafIds);
        }
        id += sizes[i];
    }

    return 0;
}

