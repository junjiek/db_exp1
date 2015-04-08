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
    IList l;
    l.insert(index);
    _list.push_back(l);
}

void Gram::sort() {
    for(auto & i : _list) {
        i.sort();
    }
}
