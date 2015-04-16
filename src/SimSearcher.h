#pragma once
#include <vector>
#include <utility>
#include <cstring>
#include <map>
#include <unordered_map>
#include <climits>
#include "Gram.h"

const int SUCCESS = 0;
const int FAILURE = 1;

#define ED 1
#define JAC 2

using namespace std;

class SimSearcher {
private:
    vector<string> _str;
    unsigned _q, _minGramSize;
    unordered_map<string, Gram> _map;
    unordered_map<string, Gram> _mapJac;
    void generateGramED(string &s, unsigned line_num);
    void generateGramJac(string &s, unsigned line_num);
public:
    SimSearcher() { _minGramSize = INT_MAX; }
    ~SimSearcher() {}
    void setQ(int q) { _q = q; }
    int createIndex(const char *filename, unsigned q);
    // calcultate T to filter the impossible ones
    int jaccardT(string &query, double threshold);
    int edT(string &query, unsigned threshold);
    //get the lists of grams for the query
    void getQueryGramList(string &query, vector<InvertedList *> &list,
                          map<int, int> &rawResult, int kind, int T);
    void scanCount(string &query, vector<InvertedList *> &list,
                   map<int, int> &rawResult, int T);
    void divideSkip(string &query, vector<InvertedList *> &list,
                    map<int, int> &rawResult, int T);
    void filter(string &query, map<int, int> &rawResult, int kind, int T);
    int levenshteinDist(string s, string t, int threshold);
    double jaccardDist(string &a, string &b, int T);
    int searchJaccard(const char *query, double threshold,
                      std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold,
                 std::vector<std::pair<unsigned, unsigned> > &result);
};

