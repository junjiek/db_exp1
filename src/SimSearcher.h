#pragma once
#include <vector>
#include <utility>
#include <cstring>
#include <map>
#include <unordered_map>
#include "Gram.h"

const int SUCCESS = 0;
const int FAILURE = 1;

#define ED 1
#define JAC 2

using namespace std;

class SimSearcher {
private:
    vector<string> _str;
    unsigned _q, _minSize;
    unordered_map<string, Gram> _map;
    void generateGram(string &s, unsigned line_num);
public:
    SimSearcher() { _minSize = 0; }
    ~SimSearcher() {}
    void setQ(int q) { _q = q; }
    int createIndex(const char *filename, unsigned q);
    // template <typename TP>
    int calcTJAC(string &query, double threshold);
    int calcTED(string &query, unsigned threshold);
    void generateList(string &query, vector<IList *> &list,
                      map<int, int> &rawResult, int kind, int T);
    void scanCount(string &query, vector<IList *> &list,
                   map<int, int> &rawResult, int T);
    void divideSkip(string &query, vector<IList *> &list,
                    map<int, int> &rawResult, int T);
    void getRawResult(string &query, map<int, int> &rawResult, int kind, int T);
    unsigned edDist(string &a, string &b, int T, unsigned threshold,
                    vector<int> &d0, vector<int> &d1);
    double jaccardDist(string &a, string &b, int T, double threshold,
                       vector<int> &d0, vector<int> &d1);
    int searchJaccard(const char *query, double threshold,
                      std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold,
                 std::vector<std::pair<unsigned, unsigned> > &result);
};

