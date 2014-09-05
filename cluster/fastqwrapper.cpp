#include "fastqwrapper.hpp"
#include <fstream>
#include <iostream>

int read_fastq(char *fname1, char *fname2, Reads &reads) {    
    std::ifstream f1(fname1),f2(fname2);
    std::string p1,p2;
    int count=0,num=0;
    while(std::getline(f1,p1) && std::getline(f2,p2)) {
        //std::cout << p1 << "\t\t\t" << p2 << "\n";
        if(count==1) {
            reads.push_back(p1);
            reads.push_back(p2);
            num++;
        }
        count+=1;
        count%=4;
    }
    return num;
}

void create_kmerhash(const Reads &reads, int k, HashTable &hashtab) {
    for (size_t j=0; j<reads.size(); ++j) {
        std::string seq = reads[j];
        int len = seq.length();
        int skipN = k;
        int64_t slid=0,mask=(1<<(2*k))-1;
        for(int i=0;i<len;i++) {
            if(seq[i]=='A') {
                slid=slid<<2;
            }
            else if(seq[i]=='C') {
                slid=slid<<2;
                slid+=1;
            }
            else if(seq[i]=='T') {
                slid=slid<<2;
                slid+=2;
            }
            else if (seq[i]=='G') {
                slid=slid<<2;
                slid+=3;
            }
            else {
                skipN=k;
                continue;
            }
            skipN--;
            if(skipN <= 0) {
                slid&=mask;
                llist* info = new llist();
                info->pos=i-k+1;
                info->readid=j;
                info->cluster = info;
                llist** point = &hashtab[slid];
                if(*point==NULL) {
                    *point=info;
                    info->next=NULL;
                }
                else {
                    info->next=*point;
                    hashtab[slid]=info;
                    //point->next=info;
                }
            }
        }
    }
    print_hashtab(hashtab);
}

void print_hashtab(const HashTable &hashtab) {
    llist* temp;
    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
        std::cout << it->first << ":";
        temp=it->second;
        while(temp!=NULL) {
            std::cout <<"--> [ " << temp->readid << "," << temp->pos <<" ] " ;
            temp=temp->next;
        }
        std::cout << "\n";
    }
}
