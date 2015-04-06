#include "SimSearcher.h"
#include <fstream>
#include <queue>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "IList.h"

#define U 0.0085

using namespace std;

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

void SimSearcher::generateGram(string &s, unsigned line_num) {
    if (s.length() < _q)
        return;
    _minSize = min(_minSize, (unsigned)s.length() - _q + 1);
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

    //sort the grams' lists
    for (auto & i : _map) {
        i.second.sort();
    }

    return (_map.empty()) ? FAILURE : SUCCESS;
}

template <typename TP>
int SimSearcher::calcT(string &query, int kind, TP threshold) {
    int T = 0;
    switch(kind) {
        case ED:
            T = max(T, (int)query.length() - (int)_q + 1 - (int)(threshold * _q));
            break;
        case JAC:
            T = max((double)threshold * ((double)query.length() - _q + 1),
                    (((double)query.length() - _q + 1) + _minSize) / (1 + 1.0 / threshold));
            break;
        default:
            break;
    }
    return T;
}

//get the lists of grams
void SimSearcher::generateList(string &query, vector<IList *> &list,
                                map<int, int> &rawResult, int kind, int T) {
    unordered_map<string, int> m;
    if (T != 0 && query.length() >= _q) {
        //get the list of q-grams
        for (int i = 0; i <= query.length() - _q; ++ i) {
            string sub = query.substr(i, _q);
            if (m.find(sub) == m.end()) {
                if (_map.find(sub) != _map.end()) {
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
        //initialize the result
        for (int i = 0; i < _str.size(); ++ i)
            rawResult[i] = 0;

        //calculate the overlap
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

void SimSearcher::divideSkip(string &query, vector<IList *> &list,
                             map<int, int> &rawResult, int T) {
    //sort the q-grams in terms of length in the descending order
    if (T == 0 || query.length() < _q)
        return;
    sort(list.begin(), list.end(), list_Compare);

    //get the L long lists
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

    //use MergeSkip to find ids that appears at leat T - L times
    while(pq.size() >= T - L) {
        vector<pair<vector<int>::iterator, vector<int>::iterator> > t;
        t.push_back(pq.top());
        pq.pop();

        while(!pq.empty()) {
            if (*(pq.top().first) == *(t[0].first)) {
                t.push_back(pq.top());
                pq.pop();
            }
            else 
                break;
        }

        int n = t.size();
        if (n >= (T - L)) {
            //check if it appears on each long list and calculate the sum of appearances
            for (auto & i : longList)
                if (i->hasKey(*(t[0].first)))
                    ++ n;

            //if it appears at least T times then add it to result
            if (n >= T)
                rawResult[*(t[0].first)] = n;

            //push next record on each popped list to the priority queue
            for (auto & i : t)
                if ((i.first + 1) != i.second)
                    pq.push(make_pair(i.first + 1, i.second));
        }
        else if (pq.size() > (T - L - 1 - n)) {
            //pop T-L-1-n smallest records from the priority queue
            for (int i = 0; i < T - L - 1 - n; ++ i) {
                t.push_back(pq.top());
                pq.pop();
            }

            //for each of the T-L-1 popped lists, find the smallest record whose value is not less than the value of the top of priority queue
            for (auto & i : t) {
                vector<int>::iterator iter =
                    lower_bound(i.first, i.second, *(pq.top().first));
                if (iter != i.second)
                    pq.push(make_pair(iter, i.second));
            }
        }
    }
}

void SimSearcher::getRawResult(string &query, map<int, int> &rawResult,
                               int kind, int T) {
    vector<IList *> list;
    rawResult.clear();
    generateList(query, list, rawResult, kind, T);
    //scanCount(query, list, rawResult, T);
    divideSkip(query, list, rawResult, T);
}

//calculate the distance
template <typename TP>
TP SimSearcher::getDistance(string &a, string &b, int T, int kind, TP threshold,
                            vector<int> &d0, vector<int> &d1) {
    double dis = 0;
    int len_a = a.length(), len_b = b.length();
    switch(kind) {
        case ED: {
            //cout << "a = " << a << endl;
            dis = threshold + 1;
            if (abs(len_a - len_b) > threshold)
                return dis;
            for (int i = 0; i <= len_a; ++ i) {
                int l = max(0, i - (int)threshold),
                            r = min(len_b, i + (int)threshold);
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
            //cout << "a = " << a << " dis = " << dis << endl;
            break;
        }
        case JAC:
            dis = 0;
            if (len_a < _q)
                return dis;
            dis = (double)T / (len_a + len_b - 2 * (_q - 1) - T);
            break;
        default:
            break;
    };
    return dis;
}

//search the similar string
template <typename TP>
int SimSearcher::searchSimilarStr(const char *query, int kind, TP threshold,
                                  vector<pair<unsigned, TP> > &result) {
    result.clear();

    //get the raw result
    string q(query);
    map<int, int> rawResult;
    getRawResult(q, rawResult, kind, calcT(q, kind, threshold));

    //eliminate the false positives
    for (auto & i : rawResult) {
        vector<int> d0(max(q.length(), _str[i.first].length()) + 1, 0);
        vector<int> d1(max(q.length(), _str[i.first].length()) + 1, 0);
        TP dis = 0;
        if (q.length() >= _str[i.first].length())
            dis = getDistance(_str[i.first], q, i.second, kind, threshold, d0, d1);
        else
            dis = getDistance(q, _str[i.first], i.second, kind, threshold, d0, d1);
        bool flag = false;
        switch(kind) {
            case ED:
                if (dis <= threshold)
                    flag = true;
                break;
            case JAC:
                if (dis >= threshold)
                    flag = true;
                break;
            default:
                break;
        }
        if (flag) {
            result.push_back(make_pair(i.first, dis));
            //cout << "fans_id = " << i.first << " " << _str[i.first] << "  " << i.second << "  " << dis << endl;
        }
    }

    return (result.empty()) ? FAILURE : SUCCESS;
}

//search the similar string in terms of Jaccard
int SimSearcher::searchJaccard(const char *query, double threshold,
                               vector<pair<unsigned, double>> &result) {
    return searchSimilarStr(query, JAC, threshold, result);
}

//search the similar string in terms of ED
int SimSearcher::searchED(const char *query, unsigned threshold,
                          vector<pair<unsigned, unsigned>> &result) {
    return searchSimilarStr(query, ED, threshold, result);
}
