#pragma once
#include <vector>
#include <utility>
#include <map>
using namespace std;

const int SUCCESS = 0;
const int FAILURE = 1;
#define BUFFSIZ 300

class SimSearcher{
private:
	int editDistance(const char* s1, const char* s2, int len1, int len2, int threshold);
	map<string, vector<pair<unsigned, unsigned>>> indexEDis;
	map<string, vector<unsigned>> indexJacc;
	// map<string,int> wordFre;
	int data_size;
	int* intersection;
	//unsigned q11;
	//vector<string> data;
	vector<int> wordNumPerID;

	void hash_init();
	double calJCD(int ind, double ths);
	unsigned calED(const char *a, int thershold, int asize,int qSiz, const char* Query);
	void createED(int lineNum, const char* s);
	void createJCD(int lineNum, const char* s);
	void defsort(int h, int t, int num);
	void mergeskip(int T, int thershold,int qSiz);
	void EDSets(int qSiz, const char* Query);
	void JCDSets(int qSiz, const char* Query);
	void print_vec_vec_int(vector<vector<int> >);

	vector<vector<int>*> data;
	vector<vector<int> > indexJCD;
	vector<vector<int> > listJCD;

	vector<int> miniStr;
	vector<const char*> dataStr;
	vector<int> lineLen;
	vector<int> visitor;
	vector<int> new_index;
	vector<int> queryCnt;

	int itemNum;
	int wordNum;
	int otherWord;
	int times;
	int minSubStrSize;
	int q;
	int querySize;
	int leave;

	int v[2][BUFFSIZ];
	int hash_v[BUFFSIZ];

public:
	SimSearcher();
	~SimSearcher();

	int createIndex(const char *filename, unsigned q);
	int searchJaccard(const char *query, double threshold, std::vector<std::pair<unsigned, double> > &result);
	int searchED(const char *query, unsigned threshold, std::vector<std::pair<unsigned, unsigned> > &result);
};

