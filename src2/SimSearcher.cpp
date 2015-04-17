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
    gramCount.clear();
    gramIdMap.clear();
    sortGramList.clear();

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
                    gramCount[gram] = 0;
                }
                // Appears > 1 times
                else {
                    num = gramCount[gram]++;
                    ostringstream sout;
                    sout << gram << num;
                    // Concat gram with num to mark apperance times of the gram 
                    // in a record
                    gram = sout.str();
                }
                if (gramIdMap.find(gram) == gramIdMap.end()) {
                    Gram newGram(gram);
                    newGram.push_back(id);
                    gramIdMap[gram] = sortGramList.size();
                    sortGramList.push_back(newGram);
                } else {
                    unsigned index = gramIdMap[gram];
                    sortGramList[index].push_back(id);
                }
            }
        }
    }
    fin.close();
    countID.resize(strings.size());
    sort(sortGramList.begin(), sortGramList.end());
    for (int i = 0; i < (int)sortGramList.size(); i++) {
        gramIdMap[sortGramList[i].getGramStr()] = i;
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
	// pair: <recordID, possible list ID>
	priority_queue<pair<unsigned, unsigned>, vector<pair<unsigned, unsigned>>,
                        heapCompare> heap;
	poppedLists.clear();
	shortResult.clear();

	unsigned topVal;
	headPos.clear();
	headPos.resize(sortGramList.size());

	// Initialize the heap
	for (int i = 0; i < shortNum; ++i)
		heap.push(make_pair(sortGramList[possibleList[i]].getList().front(), possibleList[i]));

	// MergeSkip 
	while (!heap.empty()) {
		topVal = heap.top().first;
		unsigned count = 0;
		poppedLists.clear();
        // Pop all the same value on the top
		while (!heap.empty() && heap.top().first == topVal) {
			++count;
			poppedLists.push_back(heap.top());
			heap.pop();
		}
        // Appear more than threshold times Select as a candidate
		if (count >= threshold) {
			shortResult.insert(topVal);
			countID[topVal] = count;
            vector<pair<unsigned, unsigned>>::iterator it = poppedLists.begin();
			while (it != poppedLists.end()) {
				vector<unsigned> &currList = sortGramList[it->second].getList();
				if (++headPos[it->second] < currList.size()) {
					heap.push(make_pair(currList[headPos[it->second]], it->second));
				}
                ++it;
			}
		}
        // Pop impossible ones
        else {
			for (int i = 0; !heap.empty() && i < (int)(threshold - 1 - count); ++i) {
				poppedLists.push_back(heap.top());
				heap.pop();
			}
			if (heap.empty())
				break;
			
			topVal = heap.top().first;
            vector<pair<unsigned, unsigned>>::iterator it(poppedLists.begin());
			while (it != poppedLists.end()) {
				vector<unsigned> &currList = sortGramList[it->second].getList();
				// Jump to the smallest record whose value >= heap top
				vector<unsigned>::iterator findRes = lower_bound(currList.begin(), currList.end(), topVal);
				if (findRes != currList.end()) {
					heap.push(make_pair(*findRes, it->second));
				}
				headPos[it->second] = findRes - currList.begin();
                ++it;
			}
		}
	}
}

void SimSearcher::mergeOpt(unsigned start, unsigned end, unsigned th) {
	longResult.clear();
	// MergeOpt
	for (unordered_set<unsigned>::iterator it(shortResult.begin()); it != shortResult.end(); ++it) {
        // calculate times of appearance in the L long lists.
		for (int i = start; i < (int)end; ++i) {
			// if (binary_search(sortGramList[possibleList[i]].begin(), sortGramList[possibleList[i]].end(), *it))
			if (sortGramList[possibleList[i]].hasKey(*it))
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
	unsigned num, lenQ(strlen(query)), lenW(strlen(word));
	string strQ(query), strW(word);
	for (int i = 0; i <= (int)(lenQ - q); ++i) {
		string gram(strQ.substr(i, q));
        // First appearance
		if (gramCount.find(gram) == gramCount.end()) {
			interSet.insert(gram);
			gramCount[gram] = 0;
		}
        // Appears > 1 times
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
        // First appearance
		if (gramCount.find(gram) == gramCount.end()) {
			if (interSet.find(gram) != interSet.end())
				++interNum;
			gramCount[gram] = 0;
		}
        // Appears > 1 times
		else {
			num = gramCount[gram]++;
			ostringstream sout;
			sout << gram << num;
			if (interSet.find(sout.str()) != interSet.end())
				++interNum;
			gramCount[sout.str()] = 0;
		}
	}
	
	int Gq = max(0, int(lenQ - q + 1));
    int Gw = max(0, int(lenW - q + 1));
	return double(interNum) / (Gq + Gw - interNum);
}


// Sort the id in the result
bool resultCompare(const pair<unsigned, unsigned>& a, const pair<unsigned, unsigned>& b) {
	return a.first < b.first;
}



int SimSearcher::searchJaccard(const char *query, double threshold,
                                 vector<pair<unsigned, double> > &result) {
    result.clear();
    if (divideSkip(query, jaccardT(query, threshold))) {
        // Check the candidates and 'empty'(very short) words
        double jac = 0.0;
        unordered_set<unsigned>::iterator it1 = longResult.begin();
        while (it1 != longResult.end()) {
            jac = getJac(query, strings[*it1].c_str());
            if (jac >= threshold)
                result.push_back(make_pair(*it1, jac));
            ++it1;
        }
        vector<unsigned>::iterator it2 = emptyID.begin();
        while (it2 != emptyID.end()) {
            jac = getJac(query, strings[*it2].c_str());
            if (jac >= threshold)
                result.push_back(make_pair(*it2, jac));
            ++it2;
        }
        sort(result.begin(), result.end(), resultCompare);  
    } else {
        double jac = 0.0;
        for (int i = 0; i < (int)strings.size(); ++i) {
            jac = getJac(query, strings[i].c_str());
            if (jac >= threshold)
                result.push_back(make_pair(i, jac));
        }
    }
    return (result.empty()) ? FAILURE : SUCCESS;
}

int SimSearcher::jaccardT(const char* query, double threshold) {
    int queryGramSize = max(0, (int)(strlen(query) - q + 1));
    int minGramSize = max(0, (int)(shortestStrLen - q + 1));
    return max(queryGramSize * threshold,
              (queryGramSize + minGramSize) * threshold / (1 + threshold));

}

// unsigned getED(const char *s, const char *t, int threshold) {
// 	static int distance[MAX_LEN][MAX_LEN];
// 	int slen(strlen(s)), tlen(strlen(t));
// 	 // cout << "slen" << slen << endl << tlen << endl;
// 	if (abs(slen - tlen) > threshold)
// 		return INT_MAX;

// 	for (int i = 0; i <= slen; ++i)
// 		distance[i][0] = i;
// 	for (int i = 0; i <= tlen; ++i)
// 		distance[0][i] = i;

// 	int minDist = threshold + 1;
// 	for (int i = 1; i <= slen; ++i) {
// 		int l = max(1, i - threshold);
// 		int r = min(tlen, i + threshold);
// 		minDist = threshold + 1;
// 		for (int j = l; j <= r; ++j) {
// 			if (s[i - 1] == t[j - 1])
// 				distance[i][j] = distance[i - 1][j - 1];
// 			else
// 				distance[i][j] = distance[i - 1][j - 1] + 1;
			
// 			if (abs(i - 1 - j) <= threshold && distance[i][j] > distance[i - 1][j] + 1)
// 				distance[i][j] = distance[i - 1][j] + 1;
// 			if (abs(j - 1 - i) <= threshold && distance[i][j] > distance[i][j - 1] + 1)
// 				distance[i][j] = distance[i][j - 1] + 1;
			
// 			if (distance[i][j] < minDist)
// 				minDist = distance[i][j];

// 		}
// 		if (minDist > threshold)
// 			return INT_MAX;
// 	}

//     return distance[slen][tlen];	
// }

int levenshteinDist(string s, string t, int threshold) {
    int slen = s.length();
    int tlen = t.length();

    // swap so the smaller string is t; this reduces the memory usage
    // of our buffers
    if (tlen > slen) {
        swap(s, t);
        swap(slen, tlen);
    }
    if (slen - tlen > threshold)
        return INT_MAX;

    // p is the previous and d is the current distance array; dtmp is used in swaps
    vector<int> d(tlen + 1, 0);
    vector<int> p(tlen + 1, 0);

    // the values necessary for our threshold are written; the ones after
    // must be filled with large integers since the tailing member of the threshold 
    // window in the bottom array will run min across them
    int n = 0;
    for (; n < min(tlen+1, (int)threshold + 1); ++n)
        p[n] = n;
    for (int i = n; i < tlen + 1; i++)
        p[i] = INT_MAX;
    for (int i = 0; i < tlen + 1; i++)
        d[i] = INT_MAX;

    // this is the core of the Levenshtein edit distance algorithm
    // instead of actually building the matrix, two arrays are swapped back and forth
    // the threshold limits the amount of entries that need to be computed if we're 
    // looking for a match within a set distance
    for (int row = 1; row < slen+1; ++row) {
        d[0] = row;

        // set up our threshold window
        int min_val = max(1, row - (int)threshold);
        int max_val = min(tlen+1, row + (int)threshold + 1);

        // since we're reusing arrays, we need to be sure to wipe the value left of the
        // starting index; we don't have to worry about the value above the ending index
        // as the arrays were initially filled with large integers and we progress to the right
        if (min_val > 1)
            d[min_val-1] = INT_MAX;
        int minDist = INT_MAX;
        for (int col = min_val; col < max_val; ++col) {
            if (s[row-1] == t[col-1])
                d[col] = p[col-1];
            else
                d[col] = min(p[col-1], min(d[col-1], p[col])) + 1;
            minDist = min(minDist, d[col]);
        }

        if (minDist > threshold) {
            return INT_MAX;
        }
        // swap our arrays
        swap(p, d);
    }
    return p[tlen];
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
            // ed = getED(query, strings[*it1].c_str(), threshold);
            ed = levenshteinDist(string(query), strings[*it1], threshold);

            if (ed <= threshold)
                result.push_back(make_pair(*it1, ed));
            ++it1;
        }
        vector<unsigned>::iterator it2 = emptyID.begin();
        while (it2 != emptyID.end()) {
            // ed = getED(query, strings[*it2].c_str(), threshold);
            ed = levenshteinDist(string(query), strings[*it2], threshold);
            if (ed <= threshold)
                result.push_back(make_pair(*it2, ed));
            ++it2;
        }
        sort(result.begin(), result.end(), resultCompare);
    } else {
        for (int i = 0; i < (int)strings.size(); ++i) {
            // unsigned ed = getED(query, strings[i].c_str(), threshold);
            unsigned ed = levenshteinDist(string(query), strings[i], threshold);
            if (ed <= threshold)
                result.push_back(make_pair(i, ed));
        }
    }
    return (result.empty()) ? FAILURE : SUCCESS;
}

