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
    template <typename TP>
    int calcT(string &query, int kind, TP threshold);
    void generateList(string &query, vector<IList *> &list,
                      map<int, int> &rawResult, int kind, int T);
    void scanCount(string &query, vector<IList *> &list,
                   map<int, int> &rawResult, int T);
    void divideSkip(string &query, vector<IList *> &list,
                    map<int, int> &rawResult, int T);
    void getRawResult(string &query, map<int, int> &rawResult, int kind, int T);
    template <typename TP>
    TP getDistance(string &a, string &b, int T, int kind, TP threshold,
                   vector<int> &d0, vector<int> &d1);
    template <typename TP>
    int searchSimilarStr(const char *query, int kind, TP threshold,
                         std::vector<std::pair<unsigned, TP> > &result);
    int searchJaccard(const char *query, double threshold,
                      std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold,
                 std::vector<std::pair<unsigned, unsigned> > &result);
};

