#ifndef __DSWRAPPER_HPP_
#define __DSWRAPPER_HPP_
#include <boost/pending/disjoint_sets.hpp>
#include <boost/pending/property.hpp>
class Element {
    public:
        explicit Element(int n) : mSomeInt(n) { }
        int someInt() const { return mSomeInt; }
        // The following class members are specific for the disjoint_sets
        // implementation
        size_t dsID;
        size_t dsRank;
        size_t dsParent;
    private:
        int mSomeInt;
};

inline bool operator==(Element const& lhs, Element const& rhs) {
//    std::cout << "compare : " << lhs.dsID << " vs " << rhs.dsID << "\n";
    return lhs.dsID == rhs.dsID;
}
inline bool operator!=(Element const& lhs, Element const& rhs) {
//    std::cout << "NOT ";
    return ! operator==(lhs, rhs);
}
class Parent {
    public:
        Parent(std::vector<Element>& e) : mElements(e) { }
        std::vector<Element>& mElements;
};
class MyRank {
    public:
        MyRank(std::vector<Element>& e) : mElements(e) { }
        std::vector<Element>& mElements;
};
namespace boost {
    template <>
        struct property_traits<MyRank*> {
            typedef size_t value_type;
        };
}
inline Element const& get(Parent* pa, const Element & k) {
//    std::cout << "get parent of " << k.dsID << "(" << k.dsParent << ")" << "\n";
//    std::cout << "    parent is " << pa->mElements.at(k.dsParent).dsID << "(" << pa->mElements.at(k.dsParent).dsParent << ")" << "\n";
    return pa->mElements.at(k.dsParent);
}
inline void put(Parent* pa, const Element & k, const Element & val) {
//    std::cout << "put parent " << k.dsID << "=" << val.dsID << "\n";
    pa->mElements.at(k.dsID).dsParent = val.dsID;
}
inline size_t const& get(MyRank* pa, const Element & k) {
    return pa->mElements.at(k.dsID).dsRank;
}
inline void put(MyRank* pa, Element k, size_t const& val) {
    pa->mElements.at(k.dsID).dsRank = val;
}
void printElements(std::vector<Element>& elements) {
//    std::cout << "Elements: ";
//    for (size_t i = 0; i < elements.size(); ++i) {
//        std::cout << std::setw(4) << elements[i].someInt();
//    }
//    std::cout << std::endl;
    std::cout << "ID : ";
    for (size_t i = 0; i < elements.size(); ++i) {
        std::cout << std::setw(4) << elements[i].dsID;
    }
    std::cout << std::endl;
    std::cout << "Pa : ";
    for (size_t i = 0; i < elements.size(); ++i) {
        std::cout << std::setw(4) << elements[i].dsParent;
    }
    std::cout << std::endl;
}
inline bool compareByParent(Element const& lhs, Element const& rhs) {
    return lhs.dsParent < rhs.dsParent;
}
inline bool compareBySomeInt(Element const& lhs, Element const& rhs) {
    return lhs.someInt() < rhs.someInt();
}
typedef boost::disjoint_sets<MyRank*, Parent*> DisjointSets;

#endif // __DSWRAPPER_HPP_

