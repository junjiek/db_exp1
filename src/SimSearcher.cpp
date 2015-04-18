#include "SimSearcher.h"
#include <fstream>
#include <queue>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <climits>
#include "InvertedList.h"
#define U 0.0085

using namespace std;

// clock_t start, finish;

void SimSearcher::generateGramED(string &s, unsigned line_num) {
    if (s.length() < q)
        return;
    for (int i = 0; i < s.length() - q + 1; ++ i) {
        string sub = s.substr(i, q);
        if (gramMapED.find(sub) == gramMapED.end()) {
            Gram g(sub);
            gramMapED[sub] = g;
        };
        gramMapED[sub].insert(line_num);
    }
}

void SimSearcher::generateGramJac(string &s, unsigned line_num) {
    if (s.length() == 0) return;
    unordered_set<string> subStr;
    int nend = 0;   
    int nbegin = 0;
    string sub = "";
    while(nend != -1) {
        nend = s.find(" ", nbegin);   
        if(nend == -1)
            sub = s.substr(nbegin, s.length()-nbegin);
        else  
            sub = s.substr(nbegin, nend-nbegin);
        nbegin = nend + 1;
        if (gramMapJac.find(sub) == gramMapJac.end()) {
            Gram g(sub);
            gramMapJac[sub] = g;
        };
        gramMapJac[sub].insertJac(line_num);
        subStr.insert(sub);
    }
    minSubStrSize = min(minSubStrSize, (int)subStr.size());
}

int SimSearcher::createIndex(const char *filename, unsigned q) {
    // generate q-grams
    ifstream fin(filename);
    string str;
    words.clear();
    setQ(q);
    while (getline(fin, str)) {
        generateGramED(str, words.size());
        generateGramJac(str, words.size());
        words.push_back(str);
    };
    int const initSize = 1024;
    sortedList.resize(initSize);
    candidates.resize(initSize);
    poppedLists.resize(initSize);
    headPos.resize(initSize/2);



    return (gramMapED.empty() || gramMapJac.empty()) ? FAILURE : SUCCESS;
}

void SimSearcher::getQueryGramListED(string &query) {
    unordered_map<string, int> m;
    // get the list of q-grams
    for (int i = 0; i < query.length() - q + 1; ++ i) {
        string sub = query.substr(i, q);
        if (m.find(sub) == m.end()) {
            // first-time appearance in the query grams 
            if (gramMapED.find(sub) != gramMapED.end()) {
                // has appeared in the dataset
                m[sub] = 0; 
                sortedList.push_back(&gramMapED[sub].getList(m[sub]));
            }
        }
        else {
            ++ m[sub];
            if (m[sub] < gramMapED[sub].size())
                sortedList.push_back(&gramMapED[sub].getList(m[sub]));
        }
    }
}

void SimSearcher::getQueryGramListJac(string &query) {
    for (auto & sub : querySubStr) {
        if (gramMapJac.find(sub) != gramMapJac.end()) {
            // has appeared in the dataset
            sortedList.push_back(&gramMapJac[sub].getList(0));
        }
    }
}

bool list_Compare(const InvertedList *a, const InvertedList *b) {
    return (a->size() < b->size());
};

struct heapCompare {
    bool operator() (const pair<unsigned, unsigned>& a, const pair<unsigned, unsigned>& b) {
        return a.first > b.first;
    }
};

void SimSearcher::mergeSkip(int threshold, int shortNum) {
    // pair: <recordID, possible list ID>
    priority_queue<pair<unsigned, unsigned>, vector<pair<unsigned, unsigned>>,
                        heapCompare> heap;
    poppedLists.clear();
    candidates.clear();

    int topVal;
    headPos.clear();
    headPos.resize(shortNum);

    // Initialize the heap
    for (int i = 0; i < shortNum; ++i) {
        heap.push(make_pair(sortedList[i]->getList().front(), i));
    }

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
        // Appear more than T-L times Select as a candidate
        if (count >= threshold) {
            candidates.push_back(make_pair(topVal, count));
            for (auto &pop : poppedLists) {
                vector<int> & currList = sortedList[pop.second]->getList();
                if (++headPos[pop.second] < currList.size()) {
                    heap.push(make_pair(currList[headPos[pop.second]], pop.second));
                }
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
            for (auto & pop : poppedLists) {
                vector<int> &currList = sortedList[pop.second]->getList();
                // Jump to the smallest record whose value >= heap top
                vector<int>::iterator findRes = lower_bound(currList.begin(), currList.end(), topVal);
                if (findRes != currList.end()) {
                    heap.push(make_pair(*findRes, pop.second));
                }
                headPos[pop.second] = findRes - currList.begin();
            }
        }
    }
}

void SimSearcher::mergeOpt(int begin, int end, int T) {
    // calculate times of appearance in the L long lists.
    for (auto& pair : candidates) {
        int cnt = pair.second;
        for (int i = begin; i < end; i++) {
            if (sortedList[i]->hasKey(pair.first)) {
                ++cnt;
            }
        }
        // if total appearance times >= T, add to result.
        if (cnt >= T)
            rawResult[pair.first] = true;
    }

}

void SimSearcher::divideSkip(int T) {
    //sort q-grams by length in the descending order
    sort(sortedList.begin(), sortedList.end(), list_Compare);
    //get the L longest lists
    int L = min((double(T)) / (U * log((double)(*(sortedList.back())).size()) + 1),
                double(T - 1));
    L = T - 1;
    int shortNum = sortedList.size() - int(L);
    // start = clock();
    mergeSkip(T-L, shortNum);
    // finish = clock();
    // printf("divideSkip: %.2lfs\n", (finish-start)/1000.0);
    // start = clock();

    mergeOpt(shortNum, sortedList.size(), T);
    // finish = clock();
    // printf("mergeOpt: %.2lfs\n", (finish-start)/1000.0);

    
}

void SimSearcher::filterED(string &query, int T) {
    sortedList.clear();
    rawResult.clear();
    if (T != 0 && query.length() >= q) {
        getQueryGramListED(query);
        divideSkip(T);
    } else {
        // when (T == 0 || query.length() < q) calculate directly.
        for (int i = 0; i < words.size(); ++ i)
            rawResult[i] = true;
    }
}

void SimSearcher::filterJac(string &query, int T) {
    sortedList.clear();
    rawResult.clear();
    if (T != 0 && query.length() >= q) {
        getQueryGramListJac(query);
        divideSkip(T);
    } else {
        // when (T == 0 || query.length() < q) calculate directly.
        for (int i = 0; i < words.size(); ++ i)
            rawResult[i] = true;
    }
}

double SimSearcher::jaccardDist(string& w) {
    unordered_set<string> wordSubStr;

    int nend = 0, nbegin = 0;
    while (nend != -1) {
        nend = w.find(" ", nbegin);   
        string sub = "";
        if(nend == -1)
            sub = w.substr(nbegin, w.length()-nbegin);
        else  
            sub = w.substr(nbegin, nend-nbegin);
        wordSubStr.insert(sub);
        nbegin = nend + 1;
    }

    int interNum = 0;
    for (auto & sub: querySubStr) {
        if (wordSubStr.find(sub) != wordSubStr.end())
            interNum ++;
    }

    return double(interNum) / (querySubStr.size()+wordSubStr.size()-interNum);
}

int SimSearcher::levenshteinDist(string s, string t, int threshold) {
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

int SimSearcher::jaccardT(double threshold) {
    return max(querySubStr.size() * threshold,
              (querySubStr.size() + minSubStrSize) * threshold / (1 + threshold));

}

int SimSearcher::edT(string &query, unsigned threshold) {
    return max(0, (int)query.length() - (int)q + 1 - (int)(threshold * q));
}

void SimSearcher::generateQuerySubStr(string query) {
    querySubStr.clear();
    int nend = 0, nbegin = 0;
    while (nend != -1) {
        nend = query.find(" ", nbegin);   
        string sub = "";
        if (nend == -1)
            sub = query.substr(nbegin, query.length() - nbegin);
        else
            sub = query.substr(nbegin, nend-nbegin);
       querySubStr.insert(sub);
       nbegin = nend + 1;
    }
}

//search the similar string in terms of Jaccard
int SimSearcher::searchJaccard(const char *query, double threshold,
                               vector<pair<unsigned, double>> &result) {
    result.clear();
    generateQuerySubStr(string(query));

    //get the raw result
    string queryStr(query);
    filterJac(queryStr, jaccardT(threshold));

    //eliminate false positive
    for (auto & i : rawResult) {
        double dis = jaccardDist(words[i.first]);
        if (dis >= threshold)
            result.push_back(make_pair(i.first, dis));
    }

    return SUCCESS;
}

//search the similar string in terms of ED
int SimSearcher::searchED(const char *query, unsigned threshold,
                          vector<pair<unsigned, unsigned>> &result) {
    result.clear();

    //get the raw result
    string queryStr(query);
    filterED(queryStr, edT(queryStr, threshold));

    //eliminate the false positives
    for (auto & i : rawResult) {
        unsigned dis = levenshteinDist(words[i.first], query, threshold);
        if (dis <= threshold)
            result.push_back(make_pair(i.first, dis));
    }

    return SUCCESS;
}
