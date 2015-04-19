#pragma once
#include <vector>
#include <utility>
#include <map>
#include <unordered_map>

using namespace std;

const int SUCCESS = 0;
const int FAILURE = 1;
const int T1 = 1;
const int MAXN = 1000000;
const int HASH = 95891;
const int BUFSIZE = 300;
const double EPS = 1e-8;

class SimSearcher{
private:
    int editDistance(const char* s1, const char* s2, int len1, int len2, int threshold);
    map<string, vector<pair<unsigned, unsigned>>> indexEDis;
    int* intersection;
    vector<int> wordNumPerID;

    void prepareHash();
    double calDistJac(int ind, double threshold);
    unsigned calDistED(const char *s, const char* t, int threshold);
    void createED(const char* str, int lineNum);
    void createJac(const char* str, int lineNum);
    void mysort(int b, int e, int len);
    int jaccardT(double threshold);
    int edT(unsigned threshold);
    void mergeskip(int T, int thershold);
    void getListsED(const char* query);
    void getListsJac(const char* query);
    void print_vec_vec_int(vector<vector<int> >);

    vector<vector<int>*> possibleLists;
    vector<vector<int> > wordIdxJac;
    vector<vector<int> > invertedListJac;

    vector<int> smallStr;
    vector<const char*> dataStr;
    vector<int> lineLen;
    vector<int> visitor;
    vector<int> rawResult;
    vector<int> queryIdx;

    int letterNum;
    int wordNum;
    int times;
    int minSubStrSize;
    int q;
    int querySize;
    int leave;
    int qLen;
    int v[2][BUFSIZE];
    int n_Hashq[BUFSIZE];
    unordered_map<int, vector<int>> hashED;

public:
    SimSearcher();
    ~SimSearcher();

    int createIndex(const char *filename, unsigned q);
    int searchJaccard(const char *query, double threshold, std::vector<std::pair<unsigned, double> > &result);
    int searchED(const char *query, unsigned threshold, std::vector<std::pair<unsigned, unsigned> > &result);
};

