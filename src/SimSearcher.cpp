#include "SimSearcher.h"
#include <fstream>
#include <queue>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <climits>
#include <unordered_set>
#include "InvertedList.h"
#define U 0.0085

using namespace std;

void SimSearcher::generateGramED(string &s, unsigned line_num) {
    if (s.length() < _q)
        return;
    for (int i = 0; i < s.length() - _q + 1; ++ i) {
        string sub = s.substr(i, _q);
        if (_map.find(sub) == _map.end()) {
            Gram g(sub);
            _map[sub] = g;
        };
        _map[sub].insert(line_num);
    }
}

void SimSearcher::generateGramJac(string &s, unsigned line_num) {

    int nend = 0;   
    int nbegin = 0;
    string sub = "";
    int gramSize = 0;
    while(nend != -1) {
        gramSize ++;   
        nend = s.find_first_of(" ", nbegin);   
        if(nend == -1)
            sub = s.substr(nbegin, s.length()-nbegin);
        else  
            sub = s.substr(nbegin, nend-nbegin);
        nbegin = nend + 1;
        if (_mapJac.find(sub) == _mapJac.end()) {
            Gram g(sub);
            _mapJac[sub] = g;
        };
        _mapJac[sub].insert(line_num);
    }
    _minGramSize = min(_minGramSize, gramSize);

}

int SimSearcher::createIndex(const char *filename, unsigned q) {
    // generate q-grams
    ifstream fin(filename);
    string str;
    _str.clear();
    setQ(q);
    while (getline(fin, str)) {
        generateGramED(str, _str.size());
        generateGramJac(str, _str.size());
        _str.push_back(str);
    };

    // sort the grams' lists
    for (auto & i : _map)
        i.second.sort();
    for (auto & i : _mapJac)
        i.second.sort();
    return (_map.empty() || _mapJac.empty()) ? FAILURE : SUCCESS;
}

void SimSearcher::getQueryGramList(string &query, vector<InvertedList *> &list,
                                   map<int, int> &rawResult, int kind, int T) {
    unordered_map<string, int> m;
    // get the list of q-grams
    for (int i = 0; i < query.length() - _q + 1; ++ i) {
        string sub = query.substr(i, _q);
        if (m.find(sub) == m.end()) {
            // first-time appearance in the query grams 
            if (_map.find(sub) != _map.end()) {
                // has appeared in the dataset
                m[sub] = 0;
                list.push_back(&_map[sub].getList(m[sub]));
            }
        }
        else {
            ++ m[sub];
            if (m[sub] < _map[sub].size())
                list.push_back(&_map[sub].getList(m[sub]));
        }
    }
}

void SimSearcher::scanCount(string &query, vector<InvertedList *> &list,
                            map<int, int> &rawResult, int T) {
    vector<int> counter(_str.size(), 0);

    for (auto & i : list) {
        for (auto & j : (*i).getList()) {
            counter[j] ++;
            if (counter[j] >= T)
                rawResult[j] = counter[j];
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

void SimSearcher::divideSkip(string &query, vector<InvertedList *> &list,
                             map<int, int> &rawResult, int T) {
    //sort q-grams by length in the descending order
    sort(list.begin(), list.end(), list_Compare);

    //get the L longest lists
    int L = min((double(T)) / (U * log((double)(*(list.back())).size()) + 1),
                double(T - 1));
    vector<InvertedList *> longList;
    for (int i = 0; i < L && !list.empty(); ++ i) {
        longList.push_back(list.back());
        list.pop_back();
    }

    //initialize the priority queue
    priority_queue<pair<vector<int>::iterator, vector<int>::iterator>,
                   vector<pair<vector<int>::iterator, vector<int>::iterator>>,
                   Pair_Compare> pq;
    for (auto & i : list)
        pq.push(make_pair(i->getList().begin(), i->getList().end()));

    //use MergeSkip to find ids that appear >= (T - L) times in the short lists
    while(pq.size() >= T - L) {
        vector<pair<vector<int>::iterator, vector<int>::iterator> > t;
        t.push_back(pq.top());
        pq.pop();

        while(!pq.empty()) {
            if (*(pq.top().first) == *(t[0].first)) {
                t.push_back(pq.top());
                pq.pop();
            } else 
                break;
        }

        int appearance = t.size();
        if (appearance >= (T - L)) {
            // calculate times of appearance in the L long lists.
            for (auto & i : longList)
                if (i->hasKey(*(t[0].first)))
                    ++ appearance;

            // if total appearance times >= T, add to result.
            if (appearance >= T)
                rawResult[*(t[0].first)] = appearance;

            // push next record on each popped list to the priority queue.
            for (auto & i : t)
                if ((i.first + 1) != i.second)
                    pq.push(make_pair(i.first + 1, i.second));
        } else if (pq.size() >= (T - L - appearance)) {
            // pop another (T-L-1-appearance) smallest records from the priority queue.
            for (int i = 0; i < T - L - 1 - appearance; ++ i) {
                t.push_back(pq.top());
                pq.pop();
            }

            // for the total T-L-1 popped lists, jump to the smallest record whose 
            // value >= the value of the top of pq
            for (auto & i : t) {
                vector<int>::iterator iter =
                    lower_bound(i.first, i.second, *(pq.top().first));
                if (iter != i.second)
                    pq.push(make_pair(iter, i.second));
            }
        }
    }
}

void SimSearcher::filter(string &query, map<int, int> &rawResult,
                               int kind, int T) {
    vector<InvertedList *> list;
    rawResult.clear();
    if (T != 0 && query.length() >= _q) {
        getQueryGramList(query, list, rawResult, kind, T);
        // scanCount(q uery, list, rawResult, T);
        divideSkip(query, list, rawResult, T);
    } else {
        // when (T == 0 || query.length() < _q) calculate directly.

        // initialize
        for (int i = 0; i < _str.size(); ++ i)
            rawResult[i] = 0;

        unordered_map<string, int> m;
        // calculate the overlap
        if (kind == JAC && query.length() >= _q) {
            for (int i = 0; i <= query.length() - _q; ++ i) {
                string sub = query.substr(i, _q);
                if (m.find(sub) == m.end()) {
                    if (_map.find(sub) != _map.end()) {
                        m[sub] = 0;
                        for (auto & j : _map[sub].getList(m[sub]).getList())
                            rawResult[j] ++;
                    }
                }
                else {
                    ++ m[sub];
                    if (m[sub] < _map[sub].size())
                        for (auto & j : _map[sub].getList(m[sub]).getList())
                            rawResult[j] ++;
                }
            }
        }
    }
}

double SimSearcher::jaccardDist(string& a, string& b) {
    // int len_a = a.length(), len_b = b.length();
    // if (min(len_a, len_b) < _q)
    //     return 0;
    // return (double)T / (len_a + len_b - 2 * (_q - 1) - T);
    unordered_set<string> interSet;
    int interNum = 0;

    int nend = 0, nbegin = 0;
    string sub = "";
    int gramSizeA = 0;
    while (nend != -1) {
        gramSizeA ++;   
        nend = a.find_first_of(" ", nbegin);   
        if(nend == -1)
            sub = a.substr(nbegin, a.length()-nbegin);
        else  
            sub = a.substr(nbegin, nend-nbegin);
        nbegin = nend + 1;
        interSet.insert(sub);
    }
    nend = nbegin = 0;
    int gramSizeB = 0;
    while (nend != -1) {
        gramSizeB ++;   
        nend = b.find_first_of(" ", nbegin);   
        if(nend == -1)
            sub = b.substr(nbegin, b.length()-nbegin);
        else  
            sub = b.substr(nbegin, nend-nbegin);
        nbegin = nend + 1;
        if (interSet.find(sub) != interSet.end())
            ++interNum;
    }
    return double(interNum) / (gramSizeA + gramSizeB - interNum);
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

int SimSearcher::jaccardT(string &query, double threshold) {
    double queryGramSize = (double)query.length() - _q + 1;
    return max(queryGramSize * threshold,
              (queryGramSize + _minGramSize) * threshold / (1 + threshold));

}

int SimSearcher::edT(string &query, unsigned threshold) {
    return max(0, (int)query.length() - (int)_q + 1 - (int)(threshold * _q));
}

//search the similar string in terms of Jaccard
int SimSearcher::searchJaccard(const char *query, double threshold,
                               vector<pair<unsigned, double>> &result) {
    result.clear();

    //get the raw result
    string _query(query);
    map<int, int> rawResult;
    filter(_query, rawResult, JAC, jaccardT(_query, threshold));

    //eliminate false positive
    for (auto & i : rawResult) {
        double dis = jaccardDist(_query, _str[i.first]);
        if (dis >= threshold)
            result.push_back(make_pair(i.first, dis));
    }

    return (result.empty()) ? FAILURE : SUCCESS;
}

//search the similar string in terms of ED
int SimSearcher::searchED(const char *query, unsigned threshold,
                          vector<pair<unsigned, unsigned>> &result) {
    result.clear();

    //get the raw result
    string _query(query);
    map<int, int> rawResult;
    filter(_query, rawResult, ED, edT(_query, threshold));

    //eliminate the false positives
    for (auto & i : rawResult) {
        unsigned dis = levenshteinDist(_str[i.first], _query, threshold);
        if (dis <= threshold)
            result.push_back(make_pair(i.first, dis));
    }

    return (result.empty()) ? FAILURE : SUCCESS;
}
