#ifndef GRAM_H
#define GRAM_H

#include <vector>
#include <string>
#include <iostream>
#include "IList.h"

using namespace std;

class Gram {
private:
	string _str;
	vector<IList> _list;
public:
	Gram() { _str = ""; }
	Gram(string str) { _str = str; }
	~Gram() {}
	string getString() { return _str; }
	IList & getList(int index) { return _list[index]; }
	int size() { return _list.size(); }
	void sort();
	void insert(unsigned index);
	// for debug
	void print() {
		cout << "+ " << _str << ": " << endl;
		for (auto& i : _list)
			i.print();
	}
};


#endif
