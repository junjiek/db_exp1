#include "SimSearcher.h"
#include <fstream>
#include <queue>
#include <algorithm>
#include <iostream>
#include <cassert>
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

    //sort the grams' lists
    for (auto & i : _map)
        i.second.sort();

    // for (auto & i : _map)
        // i.second.print();

    return (_map.empty()) ? FAILURE : SUCCESS;
}

void SimSearcher::getQueryGramList(string &query, vector<IList *> &list,
                                   map<int, int> &rawResult, int kind, int T) {
    unordered_map<string, int> m;
    if (T != 0 && query.length() >= _q) {
        //get the list of q-grams
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
    //sort the q-grams in terms of length in the descending order
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
    getQueryGramList(query, list, rawResult, kind, T);
    // scanCount(q uery, list, rawResult, T);
    divideSkip(query, list, rawResult, T);
}

double SimSearcher::jaccardDist(string &a, string &b, int T, double threshold,
                                vector<int> &d0, vector<int> &d1) {
    double dis = 0;
    int len_a = a.length(), len_b = b.length();

    dis = 0;
    if (len_a < _q)
        return dis;
    dis = (double)T / (len_a + len_b - 2 * (_q - 1) - T);

    return dis;
}

unsigned SimSearcher::edDist(string &a, string &b, int T, unsigned threshold,
                             vector<int> &d0, vector<int> &d1) {
    double dis = 0;
    int len_a = a.length(), len_b = b.length();
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

    return dis;
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
    getRawResult(_query, rawResult, JAC, jaccardT(_query, threshold));

    //eliminate the false positives
    for (auto & i : rawResult) {
        vector<int> d0(max(_query.length(), _str[i.first].length()) + 1, 0);
        vector<int> d1(max(_query.length(), _str[i.first].length()) + 1, 0);
        double dis = 0;
        if (_query.length() >= _str[i.first].length())
            dis = jaccardDist(_str[i.first], _query, i.second, threshold, d0, d1);
        else
            dis = jaccardDist(_query, _str[i.first], i.second, threshold, d0, d1);
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
    getRawResult(_query, rawResult, ED, edT(_query, threshold));

    //eliminate the false positives
    for (auto & i : rawResult) {
        vector<int> d0(max(_query.length(), _str[i.first].length()) + 1, 0);
        vector<int> d1(max(_query.length(), _str[i.first].length()) + 1, 0);
        unsigned dis = 0;
        if (_query.length() >= _str[i.first].length())
            dis = edDist(_str[i.first], _query, i.second, threshold, d0, d1);
        else
            dis = edDist(_query, _str[i.first], i.second, threshold, d0, d1);
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
