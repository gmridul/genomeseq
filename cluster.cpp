#include <iostream>
#include <cstdint>
#include <string>
#include <fstream>
#include <unordered_map>
using namespace std;


class llist {
    public:
    llist* next;
    char side;
    int pos;
    int entrynum;
    llist* cluster;
    llist() {
        side='c';
    }
};

int main(int argc, char*  argv[]) {
    int k,skipN=0,num=1;
    cin >> k; //k-mer size // not greater than 32
    int64_t slid,mask=(1<<(2*k))-1;
    unordered_map<int64_t,llist* > hashtab;
    ifstream f1(argv[1]),f2(argv[2]);
    string p1,p2;
    int count=0, len;

    while(getline(f1,p1) && getline(f2,p2)) {
        //f1.push_back(p1);
        //f2.push_back(p2);
        len=p1.length();
        if(count==1) {
            slid=0;
            for(int i=0;i<len-k+1;i++) {
                if(skipN>0) {
                    skipN--;
                    continue;
                }
                if(p1[i]=='A') {
                    slid=slid<<2;
                }
                else if(p1[i]=='C') {
                    slid=slid<<2;
                    slid+=1;
                }
                else if(p1[i]=='T') {
                    slid=slid<<2;
                    slid+=2;
                }
                else if (p1[i]=='G') {
                    slid=slid<<2;
                    slid+=3;
                }
                else {
                    skipN=k;
                    continue;
                }

                if(i>=k-1) {
                    slid&=mask;
                    llist* info = new llist();
                    info->side='l';
                    info->pos=i-k+1;
                    info->entrynum=num;
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

            slid=0;
            skipN=0;
            len=p2.length();
            for(int i=0;i<len-k+1;i++) {
                if(skipN>0) {
                    skipN--;
                    continue;
                }
                if(p2[i]=='A') {
                    slid=slid<<2;
                }
                else if(p2[i]=='C') {
                    slid=slid<<2;
                    slid+=1;
                }
                else if(p2[i]=='T') {
                    slid=slid<<2;
                    slid+=2;
                }
                else if(p2[i]=='G') {
                    slid=slid<<2;
                    slid+=3;
                }
                else {
                    skipN=k;
                    continue;
                }
                if(i>=k-1) {
                    slid&=mask;
                    llist* info = new llist();
                    info->side='r';
                    info->pos=i-k+1;
                    info->entrynum=num;
                    info->cluster = info;
                    llist** point = &hashtab[slid];
                    if(*point==NULL) {
                        *point=info;
                        info->next=NULL;
                    }
                    else {
                        info->next=*point;
                        hashtab[slid]=info;
                    }
                }
            }
        }

        count+=1;
        count%=4;
        num++;
    }
    // HASH TABLE CREATED
    llist* temp;
    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
        cout << it->first << ":";
        temp=it->second;
        while(temp!=NULL) {
            cout <<"--> [ " << temp->entrynum << ","<< temp->side << "," << temp->pos <<" ] " ;
            temp=temp->next;
        }
        cout << "\n";
    }

    for ( auto it = hashtab.begin(); it != hashtab.end(); ++it ) {
        temp=it->second;
        while(temp!=NULL) {

            temp=temp->next;
        }
    }
    
    



    return 0;
}
