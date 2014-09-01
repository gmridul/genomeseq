/*
 * Compile command : 
 * g++ -I/home/mridul/seqan/seqan-trunk/core/include -I/home/mridul/seqan/seqan-trunk/extras/include/ -lrt -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -std=c++0x match_reads.cpp
 */

#include <iostream>

#include "fastqwrapper.hpp"
#include "alignwrapper.hpp"
#include "dswrapper.hpp"

Reads left_reads, right_reads;

#define RANK(x) (2*x->entrynum + (x->side == 'l' ? 0 : 1) + 1) // ranks are 1, 2, ...

int match_seqs(const std::string &s1, const std::string &s2) {
    TSequence seq1 = s1, seq2 = s2;
    TAlign align;
    seqan::resize(rows(align), 2);
    seqan::assignSource(row(align, 0), seq1);
    seqan::assignSource(row(align, 1), seq2);
    int score = seqan::globalAlignment(align, seqan::Score<int,seqan::Simple>(0,-1,-1), -BOUND, BOUND);
    std::cout << "score " << score << std::endl;
    std::cout << align << std::endl;
    return score;
}

bool match_reads(llist* x,llist* y) {
    if (match_seqs(left_reads[x->entrynum], left_reads[y->entrynum]) < THRESHOLD) return false;
    if (match_seqs(right_reads[x->entrynum], right_reads[y->entrynum]) < THRESHOLD) return false;
    return true;
}

void bin_reads(const HashTable &hashtab, const std::vector<Element> &elements, DisjointSets &ds) {

    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
        std::cout << it->first << ":";
        for (auto x = it->second; (x != NULL) && (x->next != NULL); x = x->next) {
          int xrank = RANK(x);
          for (auto y = x->next; y != NULL; y = y->next) {
            int yrank = RANK(y);
            if (match_reads(x, y)) {
              ds.union_set(elements[xrank], elements[yrank]);
            }
          }
        }
    }
}

int main(int argc, char*  argv[]) {
    HashTable hashtab;
    int k = 4;
    int num = read_fastq(argv[1], argv[2], left_reads, right_reads, k, hashtab);
    print_hashtab(hashtab);

    std::vector<Element> elements;
    elements.reserve(2*num);
    for (size_t i = 0; i < elements.capacity(); ++i)
    {
        elements.push_back(Element(i+1));
    }

    for (size_t i = 0; i < elements.size(); ++i)
    {
        elements[i].dsID = i+1;
    }

    Rank rank(elements);
    Parent parent(elements);

    DisjointSets ds(&rank, &parent);

    for (size_t i = 0; i < elements.size(); ++i)
    {
        ds.make_set(elements.at(i));
    }

    bin_reads(hashtab, elements, ds);

    std::cout << "Found " << ds.count_sets(elements.begin(), elements.end()) << " sets:" << std::endl;
    printElements(elements);

    ds.compress_sets(elements.begin(), elements.end());

    std::cout << std::endl << "After path compression:" << std::endl;
    printElements(elements);

    std::sort(elements.begin(), elements.end(), compareByParent);

    std::cout << std::endl << "After path compression and sorting by parent:" << std::endl;
    printElements(elements);

    std::cout << std::endl << "Now we can iterate through all elements of each set separately using the indices:" << std::endl;
    {
        size_t first = 0;
        while (first < elements.size())
        {
            size_t currentParent = elements.at(first).dsParent;
            size_t last = first;
            while (last < elements.size() && elements.at(last).dsParent == currentParent)
            {
                ++last;
            }
            std::cout << "\tRange: [" << first << "," << last << "). Sorted elements: ";
            for (size_t i = first; i < last; ++i)
            {
                std::cout << elements.at(i).someInt() << " ";
            }
            std::cout << std::endl;
            first = last;
        }
    }
    
    return 0;
}
