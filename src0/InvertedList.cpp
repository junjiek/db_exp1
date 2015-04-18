#include "InvertedList.h"
#include <algorithm>

using namespace std;

void InvertedList::insert(unsigned index) {
    _list.push_back(index);
    _map[index] = true;
}

bool InvertedList::hasKey(int key) {
    if (_map.find(key) == _map.end())
        return false;
    return true;
}