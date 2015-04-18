#include "Gram.h"

using namespace std;

void Gram::insert(unsigned index) {
    for (auto & i : _list) {
        if (!i.hasKey(index)) {
            i.insert(index);
            return;
        }
    }
    // if the gram appeared >1 times in one string
    InvertedList l;
    l.insert(index);
    _list.push_back(l);
}

void Gram::insertJac(unsigned index) {
    if (_list.empty()) {
        InvertedList l;
        l.insert(index);
        _list.push_back(l);
        return;
    }
    if (!_list[0].hasKey(index))
        _list[0].insert(index);
}

void Gram::sort() {
    for(auto & i : _list) {
        i.sort();
    }
}
