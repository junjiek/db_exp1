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

using namespace std;

class SimSearcher {
private:
    vector<string> words;
    unsigned q;
    int minSubStrSize;
    unordered_map<string, Gram> gramMapED;
    unordered_map<string, Gram> gramMapJac;
    unordered_set<string> querySubStr;

    void setQ(int qq) { q = qq; }
    void generateGramED(string &s, unsigned line_num);
    void generateGramJac(string &s, unsigned line_num);
    void generateQuerySubStr(string query);

    // calcultate T to filter the impossible ones
    int jaccardT(double threshold);
    int edT(string &query, unsigned threshold);

    // get the lists of grams for the query
    void getQueryGramListED(string &query, vector<InvertedList *> &list);
    void getQueryGramListJac(string &query, vector<InvertedList *> &list);
    void filterED(string &query, map<int, int> &rawResult, int T);
    void filterJac(string &query, map<int, int> &rawResult, int T);
    void divideSkip(string &query, vector<InvertedList *> &list,
                    map<int, int> &rawResult, int T);

    // Calculate real dist to verify
    int levenshteinDist(string s, string t, int threshold);
    double jaccardDist(string &w);
public:
    SimSearcher() { minSubStrSize = INT_MAX; }
    ~SimSearcher() {}
    int createIndex(const char *filename, unsigned q);
    int searchJaccard(const char *query, double threshold,
                      std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold,
                 std::vector<std::pair<unsigned, unsigned> > &result);
};

