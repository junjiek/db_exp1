#include "InvertedList.h"
#include <algorithm>

using namespace std;

void InvertedList::insert(unsigned index) {
    _list.push_back(index);
    _set.insert(index);
}

bool InvertedList::hasKey(int key) {
    if (_set.find(key) == _set.end())
        return false;
    return true;
}