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
    vector<vector<int>*> possibleLists;
    vector<vector<int> > wordIdxJac;
    vector<vector<int> > invertedListJac;
    vector<int> smallStr;
    vector<const char*> dataStr;
    vector<int> lineLen;
    vector<int> visitLine;
    vector<int> rawResult;
    vector<int> queryIdx;
    unordered_map<int, vector<int>> hashED;

    int q;
    int letterNum;
    int wordNum;
    int visitor;
    int minSubStrSize;
    unsigned maxListSize;
    int querySize;
    int qLen;

    void prepareHash();
    void createED(const char* str, int lineNum);
    void createJac(const char* str, int lineNum);
    double calDistJac(int ind, double threshold);
    unsigned calDistED(const char *s, const char* t, int threshold);
    void getListsED(const char* query);
    void getListsJac(const char* query);
    int jaccardT(double threshold);
    void mysort(int b, int e, int len);
    int edT(unsigned threshold);
    void mergeskip(int T, int thershold);

public:
    SimSearcher();
    ~SimSearcher();
    int createIndex(const char *filename, unsigned q);
    int searchJaccard(const char *query, double threshold,
                      std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold,
                 std::vector<std::pair<unsigned, unsigned> > &result);
};

