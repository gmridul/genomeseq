// code to cut BAC sequence using restriction enzyme

#include <iostream>
#include <string>
#include <vector>
#include <iterator>
using namespace std;

int main() {
    string g1,g2;
    cin >> g1;
    g2.clear();
    int glen=g1.length();
    g2.resize(glen);
    for(int i=0;i<glen;i++) {
        if(g1[i]=='A') g2[i]='T';
        else if(g1[i]=='C') g2[i]='G';
        else if(g1[i]=='T') g2[i]='A';
        else if(g1[i]=='G') g2[i]='C';
    }
    string c1,c2;
    cout << g2 << endl;
    cin >> c1;
    char dir;
    cin >> dir; // l means cut left to right. r means cut right to left.
    int len=c1.length();
    size_t found = -1;
    int count=0;
    found = g1.find(c1,found+1);
    while(found!=string::npos) {
        if(dir=='r') {
            g1.insert(found+len,1,'+');
            g2.insert(found,1,'+');
            found = g1.find(c1,found+1);
            count++;
        }
        else {
            g2.insert(found+len,1,'+');
            g1.insert(found,1,'+');
            found = g1.find(c1,found+2);
            count++;
        }
    }
          
    c1=string(c1.rbegin(),c1.rend());
    glen=g1.length();
    found = g2.rfind(c1,glen-1);
    while(found!=string::npos) {
        if(dir=='r') {
            g1.insert(found+len,1,'+');
            g2.insert(found,1,'+');
            found = g2.rfind(c1,found-1);
            count++;
        }
        else {
            g2.insert(found+len,1,'+');
            g1.insert(found,1,'+');
            found = g2.rfind(c1,found-2);
            count++;
        }
    }
    
    cout << g1 << "\n" << g2 << "\n";
    return 0;
}
