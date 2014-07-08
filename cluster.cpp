#include <iostream>
#include <cstdint>
#include <string>
#include <fstream>
#include <unordered_map>
#include <list>
#include <utility> // NA NA NA NA NA BATMAN #utilitybelt
using namespace std;


class llist {
    public:
    llist* next;
    char side;
    int pos;
    llist* cluster;
    llist() {
        side='c';
    }
};

int main(int argc, char*  argv[]) {
    int k;
    cin >> k; //k-mer size
    unordered_map<char,llist* > hashtab;
    ifstream f1(argv[1]),f2(argv[2]);
    string p1,p2;
    char* slid,temp;
    unsigned long long ttk=33;
    int count=0, len;

    while(getline(f1,p1) && getline(f2,p2)) {
        if(count==1) {
            *slid=p1[0];
            temp=slid;
            for(int i=1;i<k;i++) {
                temp=temp+1;
                *temp=p1[i];
            }
            llist* info = new llist();
            info->side='l';
            info->pos=0;
            info->cluster = info;
            llist* point = hashtab[slid];
            if(point->side=='c') {
                point=info;
                info->next=NULL;
            }
            else {
                point->next=info;
                info->next=point;
            }
            len=p1.length();
            for(int i=k;i<len;i++) {
                slid.pop();
                slid.push(p1[i]);
                llist* info1 = new llist();
                info1->side='l';
                info1->pos=i;
                info1->cluster = info1;
                llist* point1 = hashtab[slid];
                if(point1->side=='c') {
                    point1=info1;
                    info1->next=NULL;
                }
                else {
                    point1->next=info1;
                    info1->next=point1;
                }
            }
            
            for(int i=0;i<k;i++) {
               slid.push(p2[i]);
            }
            llist* info2 = new llist();
            info2->side='r';
            info2->pos=0;
            info2->cluster = info2;
            llist* point2 = hashtab[slid];
            if(point2->side=='c') {
                point2=info2;
                info2->next=NULL;
            }
            else {
                point2->next=info2;
                info2->next=point2;
            }
            len=p2.length();
            for(int i=k;i<len;i++) {
                slid.pop();
                slid.push(p2[i]);
                llist* info3 = new llist();
                info3->side='l';
                info3->pos=i;
                info3->cluster = info3;
                llist* point3 = hashtab[slid];
                if(point3->side=='c') {
                    point3=info3;
                    info3->next=NULL;
                }
                else {
                    point3->next=info3;
                    info3->next=point3;
                }
            }
        }
            
        count+=1;
        count%=4;
    }
    // HASH TABLE CREATED
    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it )
        cout << " " << it->first << ":" << it->second;
    cout << endl;
    
    
    
    
    
