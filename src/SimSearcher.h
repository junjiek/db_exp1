#pragma once
#include <vector>
#include <utility>
#include <cstring>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <climits>
#include "Gram.h"

const int SUCCESS = 0;
const int FAILURE = 1;

// #define ED 1
// #define JAC 2

using namespace std;

class SimSearcher {
private:
    vector<string> _str;
    unsigned _q;
    int _minGramSize;
    unordered_map<string, Gram> _map;
    unordered_map<string, Gram> _mapJac;
    unordered_set<string> querySubStr;

    void generateGramED(string &s, unsigned line_num);
    void generateGramJac(string &s, unsigned line_num);
    void generateQuerySubStr(string query);

public:
    SimSearcher() { _minGramSize = INT_MAX; }
    ~SimSearcher() {}
    void setQ(int q) { _q = q; }
    int createIndex(const char *filename, unsigned q);
    // calcultate T to filter the impossible ones
    int jaccardT(string &query, double threshold);
    int edT(string &query, unsigned threshold);
    //get the lists of grams for the query
    void getQueryGramListED(string &query, vector<InvertedList *> &list);
    void getQueryGramListJac(string &query, vector<InvertedList *> &list);
    void divideSkip(string &query, vector<InvertedList *> &list,
                    map<int, int> &rawResult, int T);
    void filterED(string &query, map<int, int> &rawResult, int T);
    void filterJac(string &query, map<int, int> &rawResult, int T);
    int levenshteinDist(string s, string t, int threshold);
    double jaccardDist(string &w);
    int searchJaccard(const char *query, double threshold,
                      std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold,
                 std::vector<std::pair<unsigned, unsigned> > &result);
};

