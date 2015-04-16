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
#include <climits>

using namespace std;

const int MAX_LEN = 5000;
const double U = 0.0085;

SimSearcher::SimSearcher() {
	q = 0;
	shortestStrLen = MAX_LEN;
	emptyID.clear();
	strings.clear();
	sortGramList.clear();
	gramIdMap.clear();
}

// Sort from short to long.
bool gramCompare(const pair<string, vector<unsigned>>& a, const pair<string, vector<unsigned>>& b) {
	return a.second.size() < b.second.size();
}

int SimSearcher::createIndex(const char *filename, unsigned q) {
	ifstream fin(filename);
	this->q = q;

	// Create inverted-table
	unordered_map<string, vector<unsigned>> originalGram;
	vector<pair<string, vector<unsigned>>>	sortGramPair;
	gramCount.clear();
	originalGram.clear();
	sortGramPair.clear();

	string str;
	while (getline(fin, str)) {
		unsigned id = strings.size();
		strings.push_back(str);
		unsigned len = str.length();
		if (len < shortestStrLen)
			shortestStrLen = len;
		// Too short: seems empty
		if (len < q) {
			emptyID.push_back(id);
		}
		// Push back into original gram
		else {
			gramCount.clear();
			// Process same grams in one string
			unsigned num;
			for (int i = 0; i < (int)(len - q + 1); ++i) {
				string gram(str.substr(i, q));
				// First appearance
				if (gramCount.find(gram) == gramCount.end()) {
					originalGram[gram].push_back(id);
					gramCount[gram] = 0;
				}
				// Appears > 1 times
				else {
					num = gramCount[gram]++;
					ostringstream sout;
					sout << gram << num;
                    // Concat gram with num to mark apperance times of the gram 
                    // in a record
					originalGram[sout.str()].push_back(id);
					gramCount[sout.str()] = 0;
				}
			}
		}
	}
	fin.close();
	countID.resize(strings.size());

	// Sort the originalGram by the length of the id list
    unordered_map<string, vector<unsigned>>::iterator it1 = originalGram.begin();
    while (it1 != originalGram.end()) {
		sortGramPair.push_back(*it1);
        it1++;    
    }
	sort(sortGramPair.begin(), sortGramPair.end(), gramCompare);

	// Link the gram(string) with the vector index(unsigned) with unordered_map
	// and store the final sorted gram list
    vector<pair<string, vector<unsigned>>>::iterator it2 = sortGramPair.begin();
	while (it2 != sortGramPair.end()) {
		sortGramList.push_back(it2->second);
		gramIdMap[it2->first] = it2 - sortGramPair.begin();
        it2++;
	}
	maxLength = sortGramList.back().size();

	return SUCCESS;
}

void SimSearcher::getQueryGramList(const char* query) {
	// Create the possible(initial) set of gram lists
	possibleList.clear();			
	gramCount.clear();
	unsigned num  = 0;

	string queryStr(query);

	for (int i = 0; i < (int)(queryStr.length() - q + 1); ++i) {
		string gram(queryStr.substr(i, q));
        // First appearance
		if (gramCount.find(gram) == gramCount.end()) {
			unordered_map<string, unsigned>::iterator findRes;
            findRes = gramIdMap.find(gram);
			// Appeared in the input file
			if (findRes != gramIdMap.end())
				possibleList.push_back(findRes->second);
			gramCount[gram] = 0;
		}
        // Appears > 1 times
		else {
			num = gramCount[gram]++;
			ostringstream sout;
			sout << gram << num;
			unordered_map<string, unsigned>::iterator findRes;
            findRes = gramIdMap.find(sout.str());
			// Appeared in the input file
			if (findRes != gramIdMap.end())
				possibleList.push_back(findRes->second);
			gramCount[sout.str()] = 0;
		}
	}
	sort(possibleList.begin(), possibleList.end());
}

// Compare function for the heap (in MergeSkip)
struct heapCompare {
	bool operator() (const pair<unsigned, unsigned>& a, const pair<unsigned, unsigned>& b) {
		return a.first > b.first;
	}
};

void SimSearcher::mergeSkip(const char *query, unsigned threshold, int shortNum) {
	// pair: <wordID, possible gram list ID>
	priority_queue<pair<unsigned, unsigned>, vector<pair<unsigned, unsigned>>,
                        heapCompare> heap;
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
		if (cnt >= threshold) {
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
			for (int i = 0; !heap.empty() && i < (int)(threshold - 1 - cnt); ++i) {
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

void SimSearcher::mergeOpt(unsigned start, unsigned end, unsigned th) {
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

	gramCount.clear();
	/* Process same grams in one string */
	unsigned num, lenQ(strlen(query)), lenW(strlen(word));
	string strQ(query), strW(word);
	for (int i = 0; i <= (int)(lenQ - q); ++i) {
		string gram(strQ.substr(i, q));
		/* Not found: first appearance */
		if (gramCount.find(gram) == gramCount.end()) {
			interSet.insert(gram);
			gramCount[gram] = 0;
		}
		/* Not first */
		else {
			num = gramCount[gram]++;
			ostringstream sout;
			sout << gram << num;
			interSet.insert(sout.str());
			gramCount[sout.str()] = 0;
		}
	}

	gramCount.clear();
	for (int j = 0; j <= (int)(lenW - q); ++j) {
		string gram(strW.substr(j, q));
		/* Not found: first appearance */
		if (gramCount.find(gram) == gramCount.end()) {
			if (interSet.find(gram) != interSet.end())
				++interNum;
			gramCount[gram] = 0;
		}
		/* Not first */
		else {
			num = gramCount[gram]++;
			ostringstream sout;
			sout << gram << num;
			if (interSet.find(sout.str()) != interSet.end())
				++interNum;
			gramCount[sout.str()] = 0;
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
			getQueryGramList(query);

			int shortNum = possibleList.size() - int(L);

			if (shortNum > 0) {
				/* Use MergeSkip algorithm on L_short set, if not empty */	
				mergeSkip(query, T - L, shortNum);
				
				/* Use MergeOpt algorithm on L_long set. */
				mergeOpt(shortNum, possibleList.size(), L);
				
				/* Check the candidates and 'empty'(very short) words */
				double jac(0.0);
				for (unordered_set<unsigned>::iterator it(longResult.begin()); it != longResult.end(); ++it) {
					jac = getJac(query, strings[*it].c_str());
					if (jac >= threshold)
						result.push_back(make_pair(*it,jac));
				}
				for (vector<unsigned>::iterator it(emptyID.begin()); it != emptyID.end(); ++it) {
					jac = getJac(query, strings[*it].c_str());
					if (jac >= threshold)
						result.push_back(make_pair(*it, jac));
				}
				sort(result.begin(), result.end(), resultCompare);	
			} 
			else {
				double jac(0.0);
				for (int i = 0; i < (int)strings.size(); ++i) {
					jac = getJac(query, strings[i].c_str());
					if (jac >= threshold)
						result.push_back(make_pair(i, jac));
				}
			}
		}
		/* The query word is too short: Just search */
		else {
			double jac(0.0);
			for (int i = 0; i < (int)strings.size(); ++i) {
				jac = getJac(query, strings[i].c_str());
				if (jac >= threshold)
					result.push_back(make_pair(i, jac));
			}
		}
	} 
	/* Just check it one by one */
	else {
		double jac(0.0);
		for (int i = 0; i < (int)strings.size(); ++i) {
			jac = getJac(query, strings[i].c_str());
			if (jac >= threshold)
				result.push_back(make_pair(i, jac));
		}
	}

	return SUCCESS;
}

unsigned getED(const char *s, const char *t, int threshold) {
	static int distance[MAX_LEN][MAX_LEN];
	int slen(strlen(s)), tlen(strlen(t));
	 // cout << "slen" << slen << endl << tlen << endl;
	if (abs(slen - tlen) > threshold)
		return INT_MAX;

	for (int i = 0; i <= slen; ++i)
		distance[i][0] = i;
	for (int i = 0; i <= tlen; ++i)
		distance[0][i] = i;

	int minDist = threshold + 1;
	for (int i = 1; i <= slen; ++i) {
		int l = max(1, i - threshold);
		int r = min(tlen, i + threshold);
		minDist = threshold + 1;
		for (int j = l; j <= r; ++j) {
			if (s[i - 1] == t[j - 1])
				distance[i][j] = distance[i - 1][j - 1];
			else
				distance[i][j] = distance[i - 1][j - 1] + 1;
			
			if (abs(i - 1 - j) <= threshold && distance[i][j] > distance[i - 1][j] + 1)
				distance[i][j] = distance[i - 1][j] + 1;
			if (abs(j - 1 - i) <= threshold && distance[i][j] > distance[i][j - 1] + 1)
				distance[i][j] = distance[i][j - 1] + 1;
			
			if (distance[i][j] < minDist)
				minDist = distance[i][j];

		}
		if (minDist > threshold)
			return INT_MAX;
	}

    return distance[slen][tlen];	
}

int SimSearcher::jaccardT(const char* query, double threshold) {
    double queryGramSize = (double)strlen(query) - q + 1;
    return max(queryGramSize * threshold,
              (queryGramSize + shortestStrLen) * threshold / (1 + threshold));

}

int SimSearcher::edT(const char* query, unsigned threshold) {
    return strlen(query) - (int)q + 1 - (int)(threshold * q);
}

bool SimSearcher::divideSkip(const char *query, int T) {
    if (T <= 0)
        return false;
    int len = strlen(query);
    if (len <= (int)q || len < 5)
        return false;
    getQueryGramList(query);
    //get the L longest lists
    int L = min((double(T)) / (U * log((double)(maxLength)) + 1),
                double(T - 1));

    int shortNum = possibleList.size() - int(L);
    if (shortNum <= 0)
        return false;
    // Use MergeSkip on L shortest list
    mergeSkip(query, T - L, shortNum);
    mergeOpt(shortNum, possibleList.size(), T);

    return true;
}


int SimSearcher::searchED(const char *query, unsigned threshold,
                          vector<pair<unsigned, unsigned> > &result) {
    result.clear();
    if (divideSkip(query, edT(query, threshold))) {
        unsigned ed = 0;
        unordered_set<unsigned>::iterator it1 = longResult.begin();
        while (it1 != longResult.end()) {
            ed = getED(query, strings[*it1].c_str(), threshold);
            if (ed <= threshold)
                result.push_back(make_pair(*it1, ed));
            ++it1;
        }
        vector<unsigned>::iterator it2 = emptyID.begin();
        while (it2 != emptyID.end()) {
            ed = getED(query, strings[*it2].c_str(), threshold);
            if (ed <= threshold)
                result.push_back(make_pair(*it2, ed));
            ++it2;
        }
        sort(result.begin(), result.end(), resultCompare);
    } else {
        for (int i = 0; i < (int)strings.size(); ++i) {
            unsigned ed = getED(query, strings[i].c_str(), threshold);
            if (ed <= threshold)
                result.push_back(make_pair(i, ed));
        }
    }
    return (result.empty()) ? FAILURE : SUCCESS;
}

