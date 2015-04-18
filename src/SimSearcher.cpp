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

clock_t start, finish;

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

    return (gramMapED.empty() || gramMapJac.empty()) ? FAILURE : SUCCESS;
}

void SimSearcher::getQueryGramListED(string &query, vector<InvertedList *> &list) {
    unordered_map<string, int> m;
    // get the list of q-grams
    for (int i = 0; i < query.length() - q + 1; ++ i) {
        string sub = query.substr(i, q);
        if (m.find(sub) == m.end()) {
            // first-time appearance in the query grams 
            if (gramMapED.find(sub) != gramMapED.end()) {
                // has appeared in the dataset
                m[sub] = 0;
                list.push_back(&gramMapED[sub].getList(m[sub]));
            }
        }
        else {
            ++ m[sub];
            if (m[sub] < gramMapED[sub].size())
                list.push_back(&gramMapED[sub].getList(m[sub]));
        }
    }
}

void SimSearcher::getQueryGramListJac(string &query, vector<InvertedList *> &list) {
    for (auto & sub : querySubStr) {
        if (gramMapJac.find(sub) != gramMapJac.end()) {
            // has appeared in the dataset
            list.push_back(&gramMapJac[sub].getList(0));
        }
    }
}

bool list_Compare(const InvertedList *a, const InvertedList *b) {
    return (a->size() < b->size());
};

class Pair_Compare {
public:
    bool operator() (
            const pair<vector<int>::iterator, vector<int>::iterator> &a, 
            const pair<vector<int>::iterator, vector<int>::iterator> &b) {
        return *(a.first) > *(b.first);
    }
};

// void SimSearcher::mergeSkip(const char *query, unsigned threshold, int shortNum) {
//     // pair: <recordID, possible list ID>
//     priority_queue<pair<unsigned, unsigned>, vector<pair<unsigned, unsigned>>,
//                         heapCompare> heap;
//     poppedLists.clear();
//     shortResult.clear();

//     unsigned topVal;
//     headPos.clear();
//     headPos.resize(sortGramList.size());

//     // Initialize the heap
//     for (int i = 0; i < shortNum; ++i)
//         heap.push(make_pair(sortGramList[possibleList[i]].getList().front(), possibleList[i]));

//     // MergeSkip 
//     while (!heap.empty()) {
//         topVal = heap.top().first;
//         unsigned count = 0;
//         poppedLists.clear();
//         // Pop all the same value on the top
//         while (!heap.empty() && heap.top().first == topVal) {
//             ++count;
//             poppedLists.push_back(heap.top());
//             heap.pop();
//         }
//         // Appear more than threshold times Select as a candidate
//         if (count >= threshold) {
//             shortResult.insert(topVal);
//             countID[topVal] = count;
//             vector<pair<unsigned, unsigned>>::iterator it = poppedLists.begin();
//             while (it != poppedLists.end()) {
//                 vector<unsigned> &currList = sortGramList[it->second].getList();
//                 if (++headPos[it->second] < currList.size()) {
//                     heap.push(make_pair(currList[headPos[it->second]], it->second));
//                 }
//                 ++it;
//             }
//         }
//         // Pop impossible ones
//         else {
//             for (int i = 0; !heap.empty() && i < (int)(threshold - 1 - count); ++i) {
//                 poppedLists.push_back(heap.top());
//                 heap.pop();
//             }
//             if (heap.empty())
//                 break;
            
//             topVal = heap.top().first;
//             vector<pair<unsigned, unsigned>>::iterator it(poppedLists.begin());
//             while (it != poppedLists.end()) {
//                 vector<unsigned> &currList = sortGramList[it->second].getList();
//                 // Jump to the smallest record whose value >= heap top
//                 vector<unsigned>::iterator findRes = lower_bound(currList.begin(), currList.end(), topVal);
//                 if (findRes != currList.end()) {
//                     heap.push(make_pair(*findRes, it->second));
//                 }
//                 headPos[it->second] = findRes - currList.begin();
//                 ++it;
//             }
//         }
//     }
// }

void SimSearcher::mergeSkip(vector<InvertedList *> &shortLists, int shortNum,
                            vector<pair<int, int>> & candidates) {

    //initialize the priority queue
    priority_queue<pair<vector<int>::iterator, vector<int>::iterator>,
                   vector<pair<vector<int>::iterator, vector<int>::iterator>>,
                   Pair_Compare> heap;
    for (auto & i : shortLists)
        heap.push(make_pair(i->getList().begin(), i->getList().end()));
    
    //use MergeSkip to find ids that appear >= (T - L) times in the short lists
    while(heap.size() >= shortNum) {
        vector<pair<vector<int>::iterator, vector<int>::iterator> > t;
        t.push_back(heap.top());
        heap.pop();

        while(!heap.empty()) {
            if (*(heap.top().first) == *(t[0].first)) {
                t.push_back(heap.top());
                heap.pop();
            } else 
                break;
        }

        int appearance = t.size();
        if (appearance >= shortNum) {
            candidates.push_back(make_pair(*(t[0].first), appearance));
            // push next record on each popped list to the priority queue.
            for (auto & i : t)
                if ((i.first + 1) != i.second)
                    heap.push(make_pair(i.first + 1, i.second));
        } else if (heap.size() >= (shortNum - appearance)) {
            // pop another (T-L-1-appearance) smallest records from the priority queue.
            for (int i = 0; i < shortNum - 1 - appearance; ++ i) {
                t.push_back(heap.top());
                heap.pop();
            }

            // for the total T-L-1 popped lists, jump to the smallest record whose 
            // value >= the value of the top of heap
            for (auto & i : t) {
                vector<int>::iterator iter =
                    lower_bound(i.first, i.second, *(heap.top().first));
                if (iter != i.second)
                    heap.push(make_pair(iter, i.second));
            }
        }
    }
}

void SimSearcher::mergeOpt(vector<InvertedList*> &longList,
                           vector<pair<int, int>>& candidates, int T) {
    // calculate times of appearance in the L long lists.
    for (auto& pair : candidates) {
        int cnt = pair.second;
        for (auto & i : longList)
            if (i->hasKey(pair.first))
                ++ cnt;

        // if total appearance times >= T, add to result.
        if (cnt >= T)
            rawResult.insert(pair.first);
    }

}

void SimSearcher::divideSkip(vector<InvertedList *> &list, int T) {
    //sort q-grams by length in the descending order
    sort(list.begin(), list.end(), list_Compare);
    //get the L longest lists
    int L = min((double(T)) / (U * log((double)(*(list.back())).size()) + 1),
                double(T - 1));
    vector<InvertedList *> longList;
    start = clock();
    for (int i = 0; i < L && !list.empty(); ++ i) {
        longList.push_back(list.back());
        list.pop_back();
    }
    vector<pair<int, int>> candidates;
    // start = clock();
    mergeSkip(list, T-L, candidates);
    // finish = clock();
    // printf("divideSkip: %.2lfs\n", (finish-start)/1000.0);
    // start = clock();

    mergeOpt(longList, candidates, T);
    // finish = clock();
    // printf("mergeOpt: %.2lfs\n", (finish-start)/1000.0);

    
}

void SimSearcher::filterED(string &query, int T) {
    vector<InvertedList *> list;
    rawResult.clear();
    if (T != 0 && query.length() >= q) {
        getQueryGramListED(query, list);
        divideSkip(list, T);
    } else {
        // when (T == 0 || query.length() < q) calculate directly.
        for (int i = 0; i < words.size(); ++ i)
            rawResult.insert(i);
    }
}

void SimSearcher::filterJac(string &query, int T) {
    vector<InvertedList *> list;
    rawResult.clear();
    if (T != 0 && query.length() >= q) {
        getQueryGramListJac(query, list);
        divideSkip(list, T);
    } else {
        // when (T == 0 || query.length() < q) calculate directly.
        for (int i = 0; i < words.size(); ++ i)
            rawResult.insert(i);
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
        double dis = jaccardDist(words[i]);
        if (dis >= threshold)
            result.push_back(make_pair(i, dis));
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
        unsigned dis = levenshteinDist(words[i], query, threshold);
        if (dis <= threshold)
            result.push_back(make_pair(i, dis));
    }

    return SUCCESS;
}
