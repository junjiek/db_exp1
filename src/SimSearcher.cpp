#include "SimSearcher.h"
#include <fstream>
#include <queue>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <climits>
#include "IList.h"

#define U 0.0085

using namespace std;

void SimSearcher::generateGram(string &s, unsigned line_num) {
    if (s.length() < _q)
        return;
    _minGramSize = min(_minGramSize, (unsigned)s.length() - _q + 1);
    for (int i = 0; i < s.length() - _q + 1; ++ i) {
        string sub = s.substr(i, _q);
        if (_map.find(sub) == _map.end()) {
            Gram g(sub);
            _map[sub] = g;
        };
        _map[sub].insert(line_num);
    }
}

int SimSearcher::createIndex(const char *filename, unsigned q) {
    // generate q-grams
    ifstream fin(filename);
    string str;
    _str.clear();
    setQ(q);
    while (getline(fin, str)) {
        generateGram(str, _str.size());
        _str.push_back(str);
    };

    // sort the grams' lists
    for (auto & i : _map)
        i.second.sort();

    return (_map.empty()) ? FAILURE : SUCCESS;
}

void SimSearcher::getQueryGramList(string &query, vector<IList *> &list,
                                   map<int, int> &rawResult, int kind, int T) {
    unordered_map<string, int> m;
    if (T != 0 && query.length() >= _q) {
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
    else {
        // when (T == 0 || query.length() < _q) calculate directly.

        // initialize
        for (int i = 0; i < _str.size(); ++ i)
            rawResult[i] = 0;

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

void SimSearcher::scanCount(string &query, vector<IList *> &list,
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

bool list_Compare(const IList *a, const IList *b) {
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

void SimSearcher::divideSkip(string &query, vector<IList *> &list,
                             map<int, int> &rawResult, int T) {
    if (T == 0 || query.length() < _q)
        return;

    //sort q-grams by length in the descending order
    sort(list.begin(), list.end(), list_Compare);

    //get the L longest lists
    int L = min((double(T)) / (U * log((double)(*(list.back())).size()) + 1),
                double(T - 1));
    vector<IList *> longList;
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
    vector<IList *> list;
    rawResult.clear();
    getQueryGramList(query, list, rawResult, kind, T);
    // scanCount(q uery, list, rawResult, T);
    divideSkip(query, list, rawResult, T);
}

double SimSearcher::jaccardDist(string &a, string &b, int T) {
    double dis = 0;
    int len_a = a.length(), len_b = b.length();

    dis = 0;
    if (len_a < _q)
        return dis;
    dis = (double)T / (len_a + len_b - 2 * (_q - 1) - T);

    return dis;
}

unsigned SimSearcher::edDist(string &a, string &b, unsigned threshold,
                             vector<int> &d0, vector<int> &d1) {
    double dis = 0;
    int len_a = a.length(), len_b = b.length();
    //cout << "a = " << a << endl;
    dis = threshold + 1;
    if (abs(len_a - len_b) > threshold)
        return dis;
    for (int i = 0; i <= len_a; ++ i) {
        int l = max(0, i - (int)threshold);
        int r = min(len_b, i + (int)threshold);
        int minDis = threshold + 1;
        for (int j = l; j <= r; ++ j) {
            if (i == 0)
                d1[j] = j;
            else if (j == 0)
                d1[j] = i;
            else {
                if (a[i - 1] == b[j - 1])
                    d1[j] = d0[j - 1];
                else
                    d1[j] = d0[j - 1] + 1;
                if (j > l) d1[j] = min(d1[j], d1[j - 1] + 1);
                if (j < i + threshold) d1[j] = min(d1[j], d0[j] + 1);
            }
            minDis = min(minDis, d1[j]);    
        }
        if (minDis > threshold)
            return dis;
        swap(d0, d1);
    }
    dis = d0[len_b];

    return dis;
}

unsigned SimSearcher::levenshtein(string& s, string& t, unsigned threshold) {
    int slen = s.length();
    int tlen = t.length();

    // swap so the smaller string is t; this reduces the memory usage
    // of our buffers
    if(tlen > slen) {
        string stmp = s;
        s = t;
        t = stmp;
        int itmp = slen;
        slen = tlen;
        tlen = itmp;
    }

    // p is the previous and d is the current distance array; dtmp is used in swaps
    int* p = new int[tlen + 1];
    int* d = new int[tlen + 1];
    int* dtmp;

    // the values necessary for our threshold are written; the ones after
    // must be filled with large integers since the tailing member of the threshold 
    // window in the bottom array will run min across them
    int n = 0;
    for(; n < min(tlen+1, (int)threshold + 1); ++n)
        p[n] = n;
    for (int i = n; i < tlen + 1; i++)
        p[i] = INT_MAX;
    for (int i = 0; i < tlen + 1; i++)
        d[i] = INT_MAX;

    // this is the core of the Levenshtein edit distance algorithm
    // instead of actually building the matrix, two arrays are swapped back and forth
    // the threshold limits the amount of entries that need to be computed if we're 
    // looking for a match within a set distance
    for(int row = 1; row < slen+1; ++row) {
        char schar = s[row-1];
        d[0] = row;

        // set up our threshold window
        int min_val = max(1, row - (int)threshold);
        int max_val = min(tlen+1, row + (int)threshold + 1);

        // since we're reusing arrays, we need to be sure to wipe the value left of the
        // starting index; we don't have to worry about the value above the ending index
        // as the arrays were initially filled with large integers and we progress to the right
        if(min_val > 1)
            d[min_val-1] = INT_MAX;

        for(int col = min_val; col < max_val; ++col) {
            if(schar == t[col-1])
                d[col] = p[col-1];
            else 
                // min of: diagonal, left, up
                d[col] = min(p[col-1], min(d[col-1], p[col])) + 1;
        }
        // swap our arrays
        dtmp = p;
        p = d;
        d = dtmp;
    }

    if(p[tlen] == INT_MAX)
        return -1;
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
        double dis = 0;
        if (_query.length() >= _str[i.first].length())
            dis = jaccardDist(_str[i.first], _query, i.second);
        else
            dis = jaccardDist(_query, _str[i.first], i.second);
        bool flag = false;
        if (dis >= threshold)
            flag = true;
        if (flag) {
            result.push_back(make_pair(i.first, dis));
            //cout << "fans_id = " << i.first << " " << _str[i.first] << "  " << i.second << "  " << dis << endl;
        }
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
        vector<int> d0(max(_query.length(), _str[i.first].length()) + 1, 0);
        vector<int> d1(max(_query.length(), _str[i.first].length()) + 1, 0);
        unsigned dis = 0;
        cout << "-------" << endl;
        cout << _str[i.first] << endl;
        cout << _query << endl;
        cout << threshold << endl;

        if (_query.length() >= _str[i.first].length()) {
            dis = edDist(_str[i.first], _query, threshold, d0, d1);
            cout << dis << ", " << levenshtein(_str[i.first], _query, threshold) << endl;
        }
        else {
            dis = edDist(_query, _str[i.first], threshold, d0, d1);
            cout << dis << ", " << levenshtein(_query, _str[i.first], threshold) << endl;
        }
        bool flag = false;
    
        if (dis <= threshold)
            flag = true;

        if (flag) {
            result.push_back(make_pair(i.first, dis));
            //cout << "fans_id = " << i.first << " " << _str[i.first] << "  " << i.second << "  " << dis << endl;
        }
    }

    return (result.empty()) ? FAILURE : SUCCESS;
}
