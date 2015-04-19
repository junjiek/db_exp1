#pragma once
#include <vector>
#include <utility>
#include <map>
#include <unordered_map>

using namespace std;

const int SUCCESS = 0;
const int FAILURE = 1;
#define MIN(a,b,c) (a<b?(a<c?a:c):(b<c?b:c))
#define T1 1
#define MAXN 1000000
#define HASH 95891
#define BUFSIZE 300
class SimSearcher{
private:
	int editDistance(const char* s1, const char* s2, int len1, int len2, int threshold);
	map<string, vector<pair<unsigned, unsigned>>> indexEDis;
	int* intersection;
	vector<int> wordNumPerID;

	void prepareHash();
	double calJac(int ind, double ths);
	unsigned calED(const char *a, int thershold, int asize,int qLen, const char* Query);
	void createED(int lineNum, const char* str);
	void createJac(int lineNum, const char* str);
	void defsort(int h, int t, int num);
	void mergeskip(int T, int thershold,int qLen);
	void getListsED(int qLen, const char* Query);
	void getListsJac(int qLen, const char* Query);
	void print_vec_vec_int(vector<vector<int> >);

	vector<vector<int>*> rawResult;
	vector<vector<int> > wordIdxJac;
	vector<vector<int> > invertedListJac;

	vector<int> smallStr;
	vector<const char*> dataStr;
	vector<int> lineLen;
	vector<int> visitor;
	vector<int> new_index;
	vector<int> queryCnt;

	int letterNum;
	int wordNum;
	int otherWord;
	int times;
	int minSubStrSize;
	int q;
	int querySize;
	int leave;

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

