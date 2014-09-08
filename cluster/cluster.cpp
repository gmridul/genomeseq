
#include <iostream>

#include "fastqwrapper.hpp"
#include "alignwrapper.hpp"
#include "dswrapper.hpp"


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
        for (auto i = 0; i < num_leaves-1; ++i) {
            _forest.left[i] = _forest.right[i] = num_leaves+i;
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
    Clusters(int n) : num_elements(n), forest(num_elements), rank(elements), parent(elements), ds(&rank, &parent) {
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
        //assert((xp.dsParent & ~0x1) != (yp.dsParent & ~0x1));
        int nodex = set_to_node[xp.dsParent];
        int nodey = set_to_node[yp.dsParent];
        ds.link(xp, yp);
        set_to_node[ds.find_set(xp).dsParent] = forest.addNode(nodex, nodey);
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
    void cluster_reads(const Reads &reads, HashTable &hashtab) {
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
                    if ((xrank & ~0x1) == (yrank & ~0x1)) {
                        // from same pair
                        continue;
                    }
                    if(checked[xrank][yrank]) continue;
                    checked[xrank][yrank] = checked[yrank][xrank] = 1;
            
                    Element xp = ds.find_set(elements[xrank]);
                    Element yp = ds.find_set(elements[yrank]);
                    Element xp_pair = ds.find_set(elements[xrank ^ 0x1]);
                    Element yp_pair = ds.find_set(elements[yrank ^ 0x1]);
                    
                    
                    if (xp.dsID == yp.dsID || xp.dsID == yp_pair.dsID || xp_pair.dsID == yp.dsID) continue;
                    
                    //std::cout << reads[yrank] << "\n" << reverse_complement(reads[yrank]) << "\n";
                    if (match_seqs(reads[xrank], reads[yrank]) >= THRESHOLD) {
                        unite_by_parent(xp, yp);
                    }
                    else if (match_seqs(reads[xrank], reverse_complement(reads[yrank])) >= THRESHOLD) {
                        unite_by_parent(xp, yp);
                    }
                }
            }
        }
        for (int i=0; i<num_elements; i++) 
            delete [] checked[i];
        delete [] checked;
    }

    std::vector<int>& get_cluster_sizes() { return cluster_sizes; }
    uint get_root_node(int id) { return set_to_node[ds.find_set(elements[id]).dsParent]; }
    BinaryTreeInfo *getTree(int id) { return forest.getTree(set_to_node[ds.find_set(elements[id]).dsParent]); }

private:
    int num_elements;
    BinaryForest forest;
    Rank rank;
    Parent parent;
    DisjointSets ds;
    std::vector<Element> elements;
    std::vector<int> cluster_sizes;
    uint *set_to_node;
};

//define TEST_SIMPLE

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
    cls.cluster_reads(reads, hashtab);
    free_hashtab(hashtab);
#endif

    cls.finalize_clusters();
    //cls.print_elements();
    std::vector<int> sizes = cls.get_cluster_sizes();
    int id=0;
    for (size_t i = 0; i < sizes.size(); ++i) {
        std::cout << "Cluster " << i << " has " << sizes[i] << " elements\n";
        if (sizes[i] > 1) {
            //uint root =  cls.get_root_node(id);
            //std::cout << "current root = " << root << "\n";
            BinaryTreeInfo *tree = cls.getTree(id);
            //tree->print();
            //align_cluster(reads, sizes[i], tree->left, tree->right, tree->leafIds);
        } else {
            std::cout << "Singleton = " << cls.getElement(id).dsID << "\n";
        }
        id += sizes[i];
    }

    return 0;
}
