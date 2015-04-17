#ifndef GRAM
#define GRAM
#include <unordered_set>
#include <vector>
#include <string>

using namespace std;

class Gram
{
private:
    string gram;
    vector<unsigned> list;
    unordered_set<unsigned> set;
public:
    Gram(string str) { gram = str; };
    ~Gram() {};
    void push_back(unsigned index);
    bool hasKey(unsigned key);
    string getGramStr() { return gram; }
    vector<unsigned> & getList() { return list; }
    unsigned size() const { return list.size(); }
    bool operator < (const Gram& obj) const {
        return this->size() < obj.size();
    }
    bool operator > (const Gram& obj) const {
        return this->size() > obj.size();
    }
    void print();
};

#endif
