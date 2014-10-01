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

char mybwt(unsigned i, TIndex & index) { 
    typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);
    //unsigned textPos = (saAt(i, index) == 0) ? length(index) - 1 : saAt(i, index) - 1;
    //unsigned textPos = saAt(i, index);
    SAValue<TIndex>::Type p = saAt(i,index);
    if (posAtFirstLocal(p, limits)) {
        return 'B';
    } else {
        TPair p2(p.i1, p.i2-1);
        return textAt(p2, index);
    }
}

int myseqoffset(unsigned i, TIndex & index) { 
    typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);
    SAValue<TIndex>::Type p = saAt(i,index);
    return p.i2;
}

int myseqno(unsigned i, TIndex & index) { 
    typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);
    SAValue<TIndex>::Type p = saAt(i,index);
    return p.i1;
}

int main ()
{
    StringSet<CharString> myStringSet;
    //appendValue(myStringSet, "acat");
    //appendValue(myStringSet, "abd");
    appendValue(myStringSet, "alphabetagamma");
    appendValue(myStringSet, "alphadeltadelta");
    appendValue(myStringSet, "sigmaeetagamma");
    TIndex index(myStringSet);
    indexRequire(index, EsaLcp());
    indexRequire(index, EsaBwt());
    indexRequire(index, FibreSA());
    Size<TIndex>::Type minindex = 0;

    for(Size<TIndex>::Type i=0;i<length(index);++i) {
        TPair res;
        posLocalize(res, i, stringSetLimits(myStringSet));
        (i==0)?std::cout<<0:std::cout<<lcpAt(i-1,index); std::cout << " " << mybwt(i,index) << " " << res << "\n"; // getSeqNo(i,stringSetLimits(myStringSet)) << " " <<getSeqOffset(i,stringSetLimits(myStringSet)) << "\n";
    }
    for(Size<TIndex>::Type i=1;i<length(index);++i) {
        if(lcpAt(i-1,index)==((i==1)?0:lcpAt(i-2,index))) {
            if(i==length(index)-1 && minindex!=-1) {
                for(auto k=minindex;k<i-1;k++) {
                    for(auto l=k+1;l<=i-1;l++) {
                        if(mybwt(k,index)!=mybwt(l,index) || (mybwt(k,index)=='B' && mybwt(l,index)=='B')) {
                            std::cout << "< "<<myseqno(k,index)<<","<<myseqoffset(k,index)<<"> "<< "m" << "< "<<myseqno(l,index)<<","<<myseqoffset(l,index)<<"> " << "\n";
                        }
                    }
                }
            }
        }
        else if(lcpAt(i,index)<lcpAt(i-1,index)) {
            if(minindex==-1) continue;
            else {
                for(auto k=minindex;k<i-1;k++) {
                    for(auto l=k+1;l<=i-1;l++) {
                        if(mybwt(k,index)!=mybwt(l,index) || (mybwt(k,index)=='B' && mybwt(l,index)=='B')) {
                            std::cout << "< "<<myseqno(k,index)<<","<<myseqoffset(k,index)<<"> "<< suffix(indexText(index), k)<<"$" << "< "<<myseqno(l,index)<<","<<myseqoffset(l,index)<<"> " << suffix(indexText(index), l) <<"\n";
                        }
                    }
                }
                minindex=-1;
            }
        }
        else {
            minindex=i;
        }
    }
    return 0;
}

