/*
 * Compile command : 
 * g++ -I/home/mridul/seqan/seqan-trunk/core/include -I/home/mridul/seqan/seqan-trunk/extras/include/ -lrt -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -std=c++0x match_reads.cpp
 */

#include <iostream>

#include "fastqwrapper.hpp"
#include "alignwrapper.hpp"
#include "dswrapper.hpp"

Reads reads;

#define MIN_SCORE -2147483648

#define RANK(x) (x->readid) // ranks are 1, 2, ...

int match_seqs(const std::string &s1, const std::string &s2) {
    //std::cout << "aligning " << s1 << " and " << s2 << std::endl;
    std::cout.flush();
    TSequence seq1 = s1, seq2 = s2;
    TAlign align;
    seqan::resize(rows(align), 2);
    seqan::assignSource(row(align, 0), seq1);
    seqan::assignSource(row(align, 1), seq2);
    int score = seqan::globalAlignment(align, seqan::Score<int,seqan::Simple>(0,-1,-1), -BOUND, BOUND);
    if (score > MIN_SCORE) {
    std::cout << "score " << score << std::endl;
    //std::cout << align << std::endl;
    }
    return score;
}

class Clusters {
public:
    Clusters(int n) : num_elements(n), nodeid(0), rank(elements), parent(elements), ds(&rank, &parent) {
        left = new uint[num_elements];
        right = new uint[num_elements];
        set_to_node = new uint[num_elements];
        elements.reserve(num_elements);
        for (size_t i = 0; i < elements.capacity(); ++i) {
            elements.push_back(Element(i));
        }
        for (size_t i = 0; i < elements.size(); ++i) {
            elements[i].dsID = i;
            set_to_node[i] = i;
            left[i]=right[i]=num_elements+i;
        }
        for (size_t i = 0; i < elements.size(); ++i) {
            ds.make_set(elements.at(i));
        }
    }
    void unite(const Element& x, const Element& y) {
        Element xp = ds.find_set(x); 
        Element yp = ds.find_set(y); 
        if (xp != yp) {
            int nodex = set_to_node[xp.dsID];
            int nodey = set_to_node[yp.dsID];
            ds.link(xp, yp);
            left[nodeid] = nodex;
            right[nodeid] = nodey;
            set_to_node[ds.find_set(x).dsID] = num_elements + nodeid;
            nodeid++;
        }
    }
    void finalize_clusters() {
        //std::cout << "Found " << ds.count_sets(elements.begin(), elements.end()) << " sets:" << std::endl;
        //printElements(elements);

        ds.compress_sets(elements.begin(), elements.end());

        //std::cout << std::endl << "After path compression:" << std::endl;
        //printElements(elements);

        std::sort(elements.begin(), elements.end(), compareByParent);

        //std::cout << std::endl << "After path compression and sorting by parent:" << std::endl;
        //printElements(elements);

        //std::cout << std::endl << "Now we can iterate through all elements of each set separately using the indices:" << std::endl;
        size_t first = 0;
        while (first < elements.size())
        {
            size_t currentParent = elements.at(first).dsParent;
            size_t last = first;
            while (last < elements.size() && elements.at(last).dsParent == currentParent) {
                ++last;
            }
            // std::cout << "\tRange: [" << first << "," << last << "). Sorted elements: ";
            // for (size_t i = first; i < last; ++i)
            // {
            //     std::cout << elements.at(i).someInt() << " ";
            // }
            // std::cout << std::endl;
            cluster_sizes.push_back(last-first);
            first = last;
        }
    }
    void cluster_reads(HashTable &hashtab) {
        std::cout << "doing clusteing \n";
        int **checked = new int*[num_elements];
        for (int i=0; i<num_elements; i++) {
            checked[i] = new int[num_elements];
            for (int j=0; j<num_elements; j++) {
                checked[i][j]=0;
            }
        }
        for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
            std::cout << it->first << ":";
            for (auto x = it->second; (x != NULL) && (x->next != NULL); x = x->next) {
                int xrank = RANK(x);
                for (auto y = x->next; y != NULL; y = y->next) {
                    int yrank = RANK(y);
                    if (xrank > yrank) { int t=xrank; xrank=yrank; yrank=t; } 
                    //std::cout << "x=" << xrank << ", y=" << yrank << "\n";
                    if ((xrank & ~0x1) == (yrank & ~0x1)) {
                        // from same pair
                        continue;
                    }
                    if (checked[xrank][yrank]) continue;
                    if (match_seqs(reads[xrank], reads[yrank]) >= THRESHOLD) {
                        unite(elements[xrank], elements[yrank]);
                    }
                    checked[xrank][yrank] = 1;
                }
            }
        }
    }

    std::vector<Element>& get_elements() { return elements; }
    std::vector<int>& get_cluster_sizes() { return cluster_sizes; }
    uint *get_left() { return left; }
    uint *get_right() { return right; }
    uint get_root_node(int id) { return set_to_node[elements[id].dsParent]; }

private:
    int num_elements;
    Rank rank;
    Parent parent;
    DisjointSets ds;
    std::vector<Element> elements;
    std::vector<int> cluster_sizes;
    uint *left;
    uint *right;
    uint nodeid;
    uint *set_to_node;
};

char *names[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};

int main(int argc, char*  argv[]) {

    //testMuscleTree();
    //exit(0);
    //

    HashTable hashtab;
    int k = 4;
    int num = read_fastq(argv[1], argv[2], reads);
    //for (int i=0; i<reads.size(); i++) {
    //    std::cout << reads[i] << "\n";
    //}
    create_kmerhash(reads, k, hashtab);
    //print_hashtab(hashtab);

    /*
    int n = 9;
    Clusters cls(n);
    cls.unite(elements[0], elements[1]);
    cls.unite(elements[2], elements[3]);
    cls.unite(elements[2], elements[4]);
    cls.unite(elements[1], elements[4]);
    cls.unite(elements[5], elements[6]);
    cls.unite(elements[7], elements[8]);
    cls.unite(elements[5], elements[8]);
    */

    int n = 2*num;
    Clusters cls(n);
    cls.cluster_reads(hashtab);

    cls.finalize_clusters();
    std::vector<Element> elements = cls.get_elements();
    printElements(elements);
    uint *left = cls.get_left();
    uint *right = cls.get_right();
    //for (int i=0; i<n-1; ++i) {
    //    std::cout << "i=" << i << ", left[i]=" << left[i] << ", right[i]=" << right[i] << "\n";
    //}
    float *leftLength = (float *)malloc(n*sizeof(float));
    float *rightLength = (float *)malloc(n*sizeof(float));
    uint *leafIds = (uint *)malloc(n*sizeof(uint));
    char **leafNames = (char **)malloc(n*sizeof(char*));
    mseq_t *prMSeq = (mseq_t*)CKCALLOC(1,sizeof(mseq_t));
    setup_sequences(prMSeq, n, reads);
    for (int i=0; i<n; i++) {
        leftLength[i] = rightLength[i] = 1.0;
        leafIds[i] = i;
        leafNames[i] = prMSeq->sqinfo[i].name;
    }
    //testSetupSequences(prMSeq);
    std::vector<int> sizes = cls.get_cluster_sizes();
    int id=0;
    for (size_t i = 0; i < sizes.size(); ++i) {
        if (sizes[i] > 1) {
        std::cout << "Cluster " << i << " has " << sizes[i] << " elements\n";
        //  align_cluster(n, sizes[i], left, right, leftLength, rightLength, leafIds, leafNames, cls.get_root_node(id)-n, prMSeq);
        }
        id += sizes[i];
    }
    
    return 0;
}
