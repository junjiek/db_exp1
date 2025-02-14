#ifndef INVERTEDLIST_H
#define INVERTEDLIST_H

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>

using namespace std;

class InvertedList {
private:
    vector<int> _list;
    unordered_map<int, bool> _set;
public:
    InvertedList() {};
    ~InvertedList() {};
    vector<int> & getList() { return _list; }
    void sort() { std::sort(_list.begin(), _list.end()); }
    int size() const { return _list.size(); }
    void insert(unsigned index);
    bool hasKey(int key);
    // for debug
    void print() {
        cout << "  ";
        for (auto i : _list) {
            cout << i << " ";
        }
        cout << endl;
    }
};

#endif
