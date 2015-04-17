#include "Gram.h"
#include <iostream>

using namespace std;

void Gram::push_back(unsigned index) {
    if (!hasKey(index)) {
        list.push_back(index);
        set.insert(index);
    }
}

bool Gram::hasKey(unsigned key) {
    return (set.find(key) != set.end());
}

void Gram::print() {
    cout << "+ "  << gram << ": [";
    for (auto & i : list) {
        cout << i << ", ";
    }
    cout  << "]" << endl;
}