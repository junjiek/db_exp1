#pragma once
#include <vector>
#include <map>
#include <unordered_map>

using namespace std;

const int SUCCESS = 0;
const int FAILURE = 1;

const int MAXN = 1000000;
const int BUFSIZE = 300;
const int HASH = 95891;
const double EPS = 1e-8;
const double U = 0.0085;


class SimSearcher{
private:
    int q;
    int letterNum;
    int wordNum;
    int visitor;
    int minSubStrSize;
    unsigned maxListSize;
    int querySize;
    int qLen;
    unordered_map<int, vector<int>> invertedListED;
    vector<vector<int> > invertedListJac;
    vector<vector<int> > wordIdxJac;
    vector<int> visitLine;
    vector<const char*> dataStr;
    vector<int> lineLen;
    vector<int> smallStr;
    vector<vector<int>*> possibleLists;
    vector<int> rawResult;
    vector<int> queryIdxJac;

    void prepareHash();
    void buildED(const char* str, int lineNum);
    void buildJac(const char* str, int lineNum);
    double calDistJac(int ind, double threshold);
    unsigned calDistED(const char *s, const char* t, int threshold);
    void getListsED(const char* query);
    void getListsJac(const char* query);
    int jaccardT(double threshold);
    void sortListLen(int b, int e, int len);
    int edT(unsigned threshold);
    void divideSkip(int T, int thershold);

public:
    SimSearcher();
    ~SimSearcher();
    int createIndex(const char *filename, unsigned q);
    int searchJaccard(const char *query, double threshold,
                      std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold,
                 std::vector<std::pair<unsigned, unsigned> > &result);
};

