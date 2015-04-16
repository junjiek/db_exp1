#pragma once
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <queue>

const int SUCCESS = 0;
const int FAILURE = 1;

using namespace std;

class SimSearcher {
private:
    /* 'Permanent' variable */
    unsigned q;              // q: parameter of q-gram algorithm
    unsigned maxLength;      // max length of the invert table
    unsigned shortestStrLen; // shortest length of the string 

    /* vector of empty ID and sorted gram list(vector) */
    vector<unsigned> emptyID;
    vector<string>  strings;

    vector<vector<unsigned>>        sortGramList;
    unordered_map<string, unsigned> gramIdMap;

    /* 'Temporal' variable */
    vector<unsigned>                    startPos;       // start position of each list
    unordered_map<string, unsigned>     gramCount;      // handle repeated grams
    vector<unsigned>                    possibleList;   // store possible index in the sortGram.
    vector<unsigned>                    countID;        // countID[i]: appearance times of strings[i]
    unordered_set<unsigned>             shortResult;    // candidate from the 'short' part (id)
    unordered_set<unsigned>             longResult;     // candidate from the 'long' part (id)

    vector<pair<unsigned, unsigned>>    poppedLists;    // pair: <wordID, sorted gram list ID>
    // calcultate T to filter the impossible ones
    int jaccardT(const char* query, double threshold);
    int edT(const char* query, unsigned threshold);

    double getJac(const char *query, const char *word);
    void getQueryGramList(const char *query);
    void mergeSkip(const char *query, unsigned threshold, int shortNum);
    void mergeOpt(unsigned start, unsigned end, unsigned threshold);

public:
    SimSearcher() ;
    ~SimSearcher() {};

    int createIndex(const char *filename, unsigned q);
    int searchJaccard(const char *query, double threshold, std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold, std::vector<std::pair<unsigned, unsigned> > &result);
};

