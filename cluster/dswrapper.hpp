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
    return lhs.someInt() == rhs.someInt();
}
inline bool operator!=(Element const& lhs, Element const& rhs) {
    return ! operator==(lhs, rhs);
}
class Parent {
    public:
        Parent(std::vector<Element>& e) : mElements(e) { }
        std::vector<Element>& mElements;
};
class Rank {
    public:
        Rank(std::vector<Element>& e) : mElements(e) { }
        std::vector<Element>& mElements;
};
namespace boost {
    template <>
        struct property_traits<Rank*> {
            typedef size_t value_type;
        };
}
inline Element const& get(Parent* pa, Element const& k) {
    return pa->mElements.at(k.dsParent);
}
inline void put(Parent* pa, Element k, Element const& val) {
    pa->mElements.at(k.dsID).dsParent = val.dsID;
}
inline size_t const& get(Rank*, Element const& k) {
    return k.dsRank;
}
inline void put(Rank* pa, Element k, size_t const& val) {
    pa->mElements.at(k.dsID).dsRank = val;
}
void printElements(std::vector<Element>& elements) {
    std::cout << "Elements           : ";
    for (size_t i = 0; i < elements.size(); ++i) {
        std::cout << std::setw(4) << elements[i].someInt();
    }
    std::cout << std::endl;
    std::cout << "ID                 : ";
    for (size_t i = 0; i < elements.size(); ++i) {
        std::cout << std::setw(4) << elements[i].dsID;
    }
    std::cout << std::endl;
    std::cout << "Set representatives: ";
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
typedef boost::disjoint_sets<Rank*, Parent*> DisjointSets;

#endif // __DSWRAPPER_HPP_

