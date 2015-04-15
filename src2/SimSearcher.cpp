#include "SimSearcher.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <queue>
#include <cstring>
#include <cstdlib>
#include <cstdio>

using namespace std;

const int MAX_LEN = 5000;

SimSearcher::SimSearcher() {
	q = 0;
	shortestStrLen = MAX_LEN;
	emptyID.clear();
	wordList.clear();
	sortGramList.clear();
	gram2id.clear();
}

// Sort from short to long.
bool gramCompare(const pair<string, vector<unsigned>>& a, const pair<string, vector<unsigned>>& b) {
	return a.second.size() < b.second.size();
}

int SimSearcher::createIndex(const char *filename, unsigned q) {
	ifstream fin(filename);
	this->q = q;

	// Create invert-table
	unordered_map<string, vector<unsigned>> originalGram;
	vector<pair<string, vector<unsigned>>>	sortGramPair;
	countGram.clear();
	originalGram.clear();
	sortGramPair.clear();

	string str;
	while (getline(fin, str)) {
		unsigned id = wordList.size();
		wordList.push_back(str);
		unsigned len = str.length();
		if (len < shortestStrLen)
			shortestStrLen = len;
		// Too short: seems empty
		if (len < q) {
			emptyID.push_back(id);
		}
		// Push back into original gram
		else {
			countGram.clear();
			// Process same grams in one string
			unsigned num;
			for (int i = 0; i <= (int)(len - q); ++i) {
				string gram(str.substr(i, q));
				// first appearance
				if (countGram.find(gram) == countGram.end()) {
					originalGram[gram].push_back(id);
					countGram[gram] = 0;
				}
				// Not first
				else {
					num = countGram[gram]++;
					ostringstream sout;
					sout << gram << num;
					originalGram[sout.str()].push_back(id);
					countGram[sout.str()] = 0;
				}
			}
		}
	}
	fin.close();
	countID.resize(wordList.size());

	// shortestStrLen = wordList.

	// Sort the originalGram by the length of the id list(vector)
	for (unordered_map<string, vector<unsigned>>::iterator it(originalGram.begin()); it != originalGram.end(); ++it) {
		sortGramPair.push_back(*it);
	}
	sort(sortGramPair.begin(), sortGramPair.end(), gramCompare);

	// Link the gram(string) with the vector index(unsigned) with unordered_map
	// and store the final sorted gram list
	for (vector<pair<string, vector<unsigned>>>::iterator it(sortGramPair.begin()); it != sortGramPair.end(); ++it) {
		sortGramList.push_back(it->second);
		gram2id[it->first] = it - sortGramPair.begin();
	}
	maxLength = sortGramList.back().size();

	/* Check the original invert-table */
/*	ofstream logout("log.txt");
	for (unordered_map<string, vector<unsigned>>::iterator it(originalGram.begin()); it != originalGram.end(); ++it)
	{
		logout << it->first << ':';
		vector<unsigned> _vec = it->second;
		for (vector<unsigned>::iterator _it(_vec.begin()); _it != _vec.end(); ++_it)
		{
			logout << *_it << ',';
		}
		logout << endl;
	}
	logout.close();
*/
	/* Check the sorted invert-table */
/*	ofstream logsout("log_sorted.txt");
	for (vector<pair<string, vector<unsigned>>>::iterator it(sortGramPair.begin()); it != sortGramPair.end(); ++it)
	{
		logsout << (it - sortGramPair.begin()) << ' ' << it->first << ':';
		vector<unsigned> _vec = it->second;
		for (vector<unsigned>::iterator _it(_vec.begin()); _it != _vec.end(); ++_it)
		{
			logsout << *_it << ',';
		}
		logsout << endl;
	}
	logsout.close();
*/
	return SUCCESS;
}

void SimSearcher::doMakeGrams(const char* query) {
	// Create the possible(initial) set of gram lists
	possibleList.clear();			
	countGram.clear();
	unsigned num(0);

	string queryStr(query);
	int len = queryStr.length();

	for (int i = 0; i <= len - (int)q; ++i) {
		string gram(queryStr.substr(i, q));

		// Not found: first appearance
		if (countGram.find(gram) == countGram.end()) {
			unordered_map<string, unsigned>::iterator findRes = gram2id.find(gram);
			// This gram exists in the input file
			if (findRes != gram2id.end())
				possibleList.push_back(findRes->second);
			countGram[gram] = 0;
		}
		/* Not first */
		else {
			num = countGram[gram]++;
			ostringstream sout;
			sout << gram << num;
			unordered_map<string, unsigned>::iterator findRes = gram2id.find(sout.str());
			/* This gram exists in the input file */
			if (findRes != gram2id.end())
				possibleList.push_back(findRes->second);
			countGram[sout.str()] = 0;
		}
	}
	sort(possibleList.begin(), possibleList.end());
}

// Compare function for the heap (in MergeSkip algorithm )
struct heapCompare {
	bool operator() (const pair<unsigned, unsigned>& a, const pair<unsigned, unsigned>& b) {
		return a.first > b.first;
	}
};

void SimSearcher::doMergeSkip(const char *query, unsigned th, int shortNum) {
	// pair: <wordID, possible gram list ID>
	priority_queue<pair<unsigned, unsigned>, vector<pair<unsigned, unsigned>>, heapCompare> heap;
	poppedLists.clear();
	shortResult.clear();

	unsigned topVal;
	startPos.clear();
	startPos.resize(sortGramList.size());

	unsigned cnt = 0;
	/* Initialize the heap */
	for (int i = 0; i < shortNum; ++i)
		heap.push(make_pair(sortGramList[possibleList[i]].front(), possibleList[i]));

	/* MergeSkip */
	while (!heap.empty()) {
		topVal = heap.top().first;
		cnt = 0;
		poppedLists.clear();
		while (!heap.empty() && heap.top().first == topVal) {
			++cnt;
			poppedLists.push_back(heap.top());
			heap.pop();
		}
		if (cnt >= th) {
			shortResult.insert(topVal);
			countID[topVal] = cnt;
			for (vector<pair<unsigned, unsigned>>::iterator it(poppedLists.begin()); it != poppedLists.end(); ++it) {
				vector<unsigned> &currList = sortGramList[it->second];
				if (++startPos[it->second] < currList.size()) {
					heap.push(make_pair(currList[startPos[it->second]], it->second));
				}
			}
		}
		else {
			for (int i = 0; !heap.empty() && i < (int)(th - 1 - cnt); ++i) {
				poppedLists.push_back(heap.top());
				heap.pop();
			}
			if (heap.empty())
				break;
			
			topVal = heap.top().first;
			for (vector<pair<unsigned, unsigned>>::iterator it(poppedLists.begin()); it != poppedLists.end(); ++it) {
				vector<unsigned> &currList = sortGramList[it->second];
				
				vector<unsigned>::iterator findRes = lower_bound(currList.begin(), currList.end(), topVal);
				if (findRes != currList.end()) {
					heap.push(make_pair(*findRes, it->second));
				}
				startPos[it->second] = findRes - currList.begin();
				/*bool brk(0);
				for (vector<unsigned>::iterator jt(currList.begin()); jt != currList.end(); ++jt)
				{
					if (*jt >= topVal)
					{
						heap.push(make_pair(*jt, it->second));
						// startPos[it->second] = jt - currList.begin();
						startPos[it->second] = jt - currList.begin();
						brk = 1;
						break;
					}
				}
				if (!brk)
				{
					startPos[it->second] =currList.size();
				}*/
			}
		}
	}

/*	cout << "start pos: " << endl;
	for (int i = 0; i < shortNum; ++i)
	{
		cout << i << ',' << startPos[i] << endl;
	}
*/
}

void SimSearcher::doMergeOpt(unsigned start, unsigned end, unsigned th) {
	longResult.clear();
	/* MergeOpt */
	for (unordered_set<unsigned>::iterator it(shortResult.begin()); it != shortResult.end(); ++it) {
		/* For each 'long' lists */
		for (int i = start; i < (int)end; ++i) {
			if (binary_search(sortGramList[possibleList[i]].begin(), sortGramList[possibleList[i]].end(), *it))
				++countID[*it];
		}
		if (countID[*it] >= th)
			longResult.insert(*it);
	}
}

double SimSearcher::getJac(const char *query, const char *word) {
	unordered_set<string> interSet;
	int interNum(0);

	countGram.clear();
	/* Process same grams in one string */
	unsigned num, lenQ(strlen(query)), lenW(strlen(word));
	string strQ(query), strW(word);
	for (int i = 0; i <= (int)(lenQ - q); ++i) {
		string gram(strQ.substr(i, q));
		/* Not found: first appearance */
		if (countGram.find(gram) == countGram.end()) {
			interSet.insert(gram);
			countGram[gram] = 0;
		}
		/* Not first */
		else {
			num = countGram[gram]++;
			ostringstream sout;
			sout << gram << num;
			interSet.insert(sout.str());
			countGram[sout.str()] = 0;
		}
	}

	countGram.clear();
	for (int j = 0; j <= (int)(lenW - q); ++j) {
		string gram(strW.substr(j, q));
		/* Not found: first appearance */
		if (countGram.find(gram) == countGram.end()) {
			if (interSet.find(gram) != interSet.end())
				++interNum;
			countGram[gram] = 0;
		}
		/* Not first */
		else {
			num = countGram[gram]++;
			ostringstream sout;
			sout << gram << num;
			if (interSet.find(sout.str()) != interSet.end())
				++interNum;
			countGram[sout.str()] = 0;
		}
	}
	
	int Gq(max(0, int(lenQ - q + 1))), Gw(max(0, int(lenW - q + 1)));
	return double(interNum) / (Gq + Gw - interNum);
}


// Sort the id in the result
bool resultCompare(const pair<unsigned, unsigned>& a, const pair<unsigned, unsigned>& b) {
	return a.first < b.first;
}


int SimSearcher::searchJaccard(const char *query, double threshold, vector<pair<unsigned, double> > &result) {
	result.clear();

	int Gq(max(0, int(strlen(query) - q + 1))), Gs(max(0, int(shortestStrLen - q + 1)));
	int T = (int)max(threshold * Gq, (Gq + Gs) * threshold / (1 + threshold));

	const double mu = 0.0085;

	// T  = -1;
	// cout << "T = " << T << endl;

	/* Using DivideSkip algorithm */
	if (T > 0) {
		unsigned L = T / (mu * log10(double(maxLength)) + 1);		// important parameter in the DivideSkip algorithm

		// cout << "L = " << L << endl;
		unsigned len = strlen(query);
		/* Parse the grams if the query string is long enough */
		if (len > q && len >= 10) {
			doMakeGrams(query);

			int shortNum = possibleList.size() - int(L);

			if (shortNum > 0) {
				/* Use MergeSkip algorithm on L_short set, if not empty */	
				doMergeSkip(query, T - L, shortNum);
				
				/* Use MergeOpt algorithm on L_long set. */
				doMergeOpt(shortNum, possibleList.size(), L);
				
				/* Check the candidates and 'empty'(very short) words */
				double jac(0.0);
				for (unordered_set<unsigned>::iterator it(longResult.begin()); it != longResult.end(); ++it) {
					jac = getJac(query, wordList[*it].c_str());
					if (jac >= threshold)
						result.push_back(make_pair(*it,jac));
				}
				for (vector<unsigned>::iterator it(emptyID.begin()); it != emptyID.end(); ++it) {
					jac = getJac(query, wordList[*it].c_str());
					if (jac >= threshold)
						result.push_back(make_pair(*it, jac));
				}
				sort(result.begin(), result.end(), resultCompare);	
			} 
			else {
				double jac(0.0);
				for (int i = 0; i < (int)wordList.size(); ++i) {
					jac = getJac(query, wordList[i].c_str());
					if (jac >= threshold)
						result.push_back(make_pair(i, jac));
				}
			}
		}
		/* The query word is too short: Just search */
		else {
			double jac(0.0);
			for (int i = 0; i < (int)wordList.size(); ++i) {
				jac = getJac(query, wordList[i].c_str());
				if (jac >= threshold)
					result.push_back(make_pair(i, jac));
			}
		}
	} 
	/* Just check it one by one */
	else {
		double jac(0.0);
		for (int i = 0; i < (int)wordList.size(); ++i) {
			jac = getJac(query, wordList[i].c_str());
			if (jac >= threshold)
				result.push_back(make_pair(i, jac));
		}
	}

	return SUCCESS;
}

unsigned getED(const char *query, const char *word, int th) {
	static int distance[MAX_LEN][MAX_LEN];

	int lenQ(strlen(query)), lenW(strlen(word));
	 // cout << "lenQ" << lenQ << endl << lenW << endl;
	int ed(0);
	if (abs(lenQ - lenW) > th)
		return th + 1;

	for (int q = 0; q <= lenQ; ++q)
		distance[q][0] = q;
	for (int w = 0; w <= lenW; ++w)
		distance[0][w] = w;

	int wSt, wEd, tmpMin(th + 1);
	for (int q = 1; q <= lenQ; ++q) {
		wSt = max(1, q - th);
		wEd = min(lenW, q + th);
		tmpMin = th + 1;
		for (int w = wSt; w <= wEd; ++w) {
			// query[q - 1] is the q_th char of query
			if (query[q - 1] == word[w - 1])
				distance[q][w] = distance[q - 1][w - 1];
			else
				distance[q][w] = distance[q - 1][w - 1] + 1;
			
			// distance[q][w] = distance[q - 1][w - 1] + (query[q - 1] != word[w - 1]?1:0);
			if (abs(q - 1 - w) <= th && distance[q][w] > distance[q - 1][w] + 1)
				distance[q][w] = distance[q - 1][w] + 1;
			if (abs(w - 1 - q) <= th && distance[q][w] > distance[q][w - 1] + 1)
				distance[q][w] = distance[q][w - 1] + 1;
			
			if (distance[q][w] < tmpMin)
				tmpMin = distance[q][w];

			// cout << q << ',' << w << ':' << distance[q][w] << ' ';
		}
		// cout << endl;
		// Break: no need to search any more
		if (tmpMin > th)
			break;
	}
	if (tmpMin <= th)
		ed = distance[lenQ][lenW];
	else
		ed = th + 1;
	
	return ed;
}

int SimSearcher::searchED(const char *query, unsigned threshold, vector<pair<unsigned, unsigned> > &result) {
	result.clear();
	/* T-occurrence threshold */
	int T = strlen(query) - q + 1 - threshold * q;
	const double mu = 0.0085;

	// T  = -1;
	// cout << "T = " << T;

	/* Using DivideSkip algorithm */
	if (T > 0) {
		unsigned L = T / (mu * log10(double(maxLength)) + 1);		// important parameter in the DivideSkip algorithm

		// cout << "; L = " << L << endl;
		unsigned len = strlen(query);
		/* Parse the grams if the query string is long enough */
		if (len > q && len >= 5) {
			doMakeGrams(query);

			int shortNum = possibleList.size() - int(L);

			if (shortNum > 0) {
				/* Use MergeSkip algorithm on L_short set, if not empty */	
				doMergeSkip(query, T - L, shortNum);
				
				/* Use MergeOpt algorithm on L_long set. */
				doMergeOpt(shortNum, possibleList.size(), T);
				
				/* Check the candidates and 'empty'(very short) words */
				unsigned ed = 0;
				for (unordered_set<unsigned>::iterator it(longResult.begin()); it != longResult.end(); ++it) {
					ed = getED(query, wordList[*it].c_str(), threshold);
					if (ed <= threshold)
						result.push_back(make_pair(*it, ed));
				}
				for (vector<unsigned>::iterator it(emptyID.begin()); it != emptyID.end(); ++it) {
					ed = getED(query, wordList[*it].c_str(), threshold);
					if (ed <= threshold) {
						result.push_back(make_pair(*it, ed));
					}
				}
				sort(result.begin(), result.end(), resultCompare);
			} 
			else {
				for (int i = 0; i < (int)wordList.size(); ++i) {
					unsigned ed = getED(query, wordList[i].c_str(), threshold);
					if (ed <= threshold)
						result.push_back(make_pair(i, ed));
				}
			}
		}
		/* The query word is too short: Just search */
		else {
			for (int i = 0; i < (int)wordList.size(); ++i) {
				unsigned ed = getED(query, wordList[i].c_str(), threshold);
				if (ed <= threshold)
					result.push_back(make_pair(i, ed));
			}
		}
	} 
	/* Just check it one by one */
	else {
		for (int i = 0; i < (int)wordList.size(); ++i) {
			unsigned ed = getED(query, wordList[i].c_str(), threshold);
			if (ed <= threshold)
				result.push_back(make_pair(i, ed));
		}
	}

	return SUCCESS;
}

