#include "SimSearcher.h"

#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>

#include <unordered_set>
#include <climits>

using namespace std;

int wordExis[MAXN][BUFSIZE];
int globalWordIdx[MAXN];
int visit[MAXN];

SimSearcher::SimSearcher() {
    minSubStrSize = INT_MAX;
    wordNum = 0;
    memset(globalWordIdx, -1, sizeof(globalWordIdx));
    letterNum = 0;
    dataStr.clear();
    hashED.clear();
    smallStr.clear();
}

SimSearcher::~SimSearcher() {}

void SimSearcher::mysort(int b, int e, int len) {
    int i = b, j = e;
    unsigned pivot = rawResult[(b + e)/2]->size();
    do {
        while (rawResult[i]->size() > pivot)
            i++;
        while (rawResult[j]->size() < pivot)
            j--;
        if (i <= j) {
            if (rawResult[i]->size() != rawResult[j]->size())
                swap(rawResult[i], rawResult[j]);
            i++;
            j--;
        }
    } while (i <= j);
    if (j > len)
        mysort(b, j, len);
    if (i <= len)
        mysort(i, e, len);
}

void SimSearcher::mergeskip(int T, int thershold) {
    if (T < 1) {
        int j = dataStr.size() - 1;
        for (int i = smallStr.size()-1; i >= 0; --i) {
            for (int k = j; k > smallStr[i]; --k)
                if (abs(lineLen[k]-qLen) <= thershold)
                    newIdx.push_back(k);
            j = smallStr[i] - 1;
        }
        for (int k = j; k >= 0; --k)
            if (abs(lineLen[k]-qLen) <= thershold)
                newIdx.push_back(k);
        return;
    }
    ++times;
    int occur = T1;
    int len = rawResult.size();
    leave = T - occur;
    if (leave < len && leave > 0)
        mysort(0, len-1, leave - 1);
    int i = leave;
    while(i < len) {
        vector<int> &curr = *(rawResult[i]);
        for (int j = curr.size() - 1; j >= 0; --j) {
            int temp = curr[j];
            if (visitor[temp] != times) {
                visitor[temp] = times;
                if (abs(lineLen[temp] - qLen) <= thershold)
                    newIdx.push_back(temp);
            }
        }
        i++;
    }
}
double SimSearcher::calDistJac(int index, double thershold) {
    vector<int> &wordIdx = wordIdxJac[index];
    int length = wordIdx.size(), bsize = querySize;
    if (length * thershold > bsize || bsize * thershold > length)
        return 0;
    int intersec = 0, q = otherWord + length;
    int i = 0, j = 0;
    while (i < length) {
        while (wordIdx[i] > queryCnt[j]) {
            ++j;
            ++q;
        }
        if (wordIdx[i++] == queryCnt[j]) {
            ++intersec;
            ++j;
        }
    }
    return (double)intersec / (q + bsize - j);
}

unsigned SimSearcher::calDistED(const char *s, const char *t, int threshold) {
    static int distance[BUFSIZE][BUFSIZE];
    int slen(strlen(s)), tlen(strlen(t));
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


void SimSearcher::createED(const char * str, int lineNum) {
    if (lineLen[lineNum] < q) {
        smallStr.push_back(lineNum);
        return;
    }
    int hashCode = 0;
    for (int i = 0; i < q; i++) {
        hashCode = hashCode * HASH + str[i];
    }
    hashED[hashCode].push_back(lineNum);
    for (int i = q; i < lineLen[lineNum]; i++) {
        hashCode = hashCode * HASH - n_Hashq[(int)(str[i-q])] + str[i];
        vector<int> &list = hashED[hashCode];
        if (list.empty() || list.back() != lineNum) {
            list.push_back(lineNum);
        }
    }
}

void SimSearcher::createJac(const char * str, int lineNum) {
    vector<int> wordIdx;
    wordIdx.clear();
    int curr = 0;
    int subStrSize = 0;

    for (int i = 0 ; i <= lineLen[lineNum]; i++) {
        if (i == lineLen[lineNum] || str[i] == ' ') {
            // The word has never appeared in dataset before
            if (globalWordIdx[curr] == -1) {
                globalWordIdx[curr] = wordNum;
                wordNum ++;
                vector<int> newWordList;
                invertedListJac.push_back(newWordList);
            }
            // The word appear first time in str
            int idx = globalWordIdx[curr];
            if (invertedListJac[idx].empty() || invertedListJac[idx].back() != lineNum) {
                invertedListJac[idx].push_back(lineNum);
                wordIdx.push_back(idx);
                subStrSize ++;
            }
            curr = 0;

        } else {
            int &next = wordExis[curr][(int)str[i]];
            if (next == 0) {
                next = ++letterNum;
            }
            curr = next;
        }
    }
    if (minSubStrSize > subStrSize) {
        minSubStrSize = subStrSize;
    }
    sort(wordIdx.begin(), wordIdx.end());
    wordIdxJac.push_back(wordIdx);
}

inline int mypow(int x, int y) {
    int result = 1;
    for (int i = 0; i < y; i ++)
        result *= x;
    return result;
}

void SimSearcher::prepareHash() {
    int Hashq = mypow(HASH, q);
    for (int n = 0; n < BUFSIZE; n++) {
        n_Hashq[n] = n * Hashq;
    }
}

int SimSearcher::createIndex(const char *filename, unsigned q) {
    this->q = q;
    prepareHash();
    ifstream fin(filename);
    string line;
    char * buf;
    while (getline(fin, line)) {
        lineLen.push_back((int)line.length());
        int lineNum = dataStr.size();
        buf = (char*)malloc(1000);
        strcpy(buf, line.c_str()); 
        dataStr.push_back(buf);
        createED(buf, lineNum);
        createJac(buf, lineNum);
    };
    fin.close();
    visitor.resize(dataStr.size());
    return SUCCESS;

}

void SimSearcher::getListsED(const char* query) {
    rawResult.clear();
    querySize = qLen + 1 - q;
    if (qLen < q) return;
    int hashCode = 0;
    for (int i = 0; i < q; i++) {
        hashCode = hashCode * HASH + query[i];
    }
    unordered_map<int, vector<int>>::iterator iter = hashED.find(hashCode);
    if (iter != hashED.end()) {
        rawResult.push_back(&(iter->second));
    }
    for (int i = q; i < qLen; i++) {
        hashCode = hashCode * HASH - n_Hashq[(int)(query[i-q])] + query[i];
        iter = hashED.find(hashCode);
        if (iter != hashED.end()) {
            rawResult.push_back(&(iter->second));
        }
    }
}
void SimSearcher::getListsJac(const char* query) {
    rawResult.clear();
    otherWord = 0;
    queryCnt.clear();
    int curr = 0;
    ++times;
    bool find = false;
    for (int i = 0; i <= qLen; ++i) {
        if (i == qLen || query[i] == ' ') {
            int num = globalWordIdx[curr];
            if (!find && num != -1) {
                rawResult.push_back(&invertedListJac[num]);
                if (visit[num] != times) {
                    visit[num] = times;
                    queryCnt.push_back(num);
                }
                //cout << "push" << idx << endl;
            }
            else
                ++otherWord;
            curr = 0;
            find = false;
        } else {
            if (find)
                continue;
            int to = query[i];
            int &next = wordExis[curr][to];
            if (next == 0) {
                find = true;
                continue;
            }
            curr = next;
        }
    }
    sort(queryCnt.begin(), queryCnt.end());
    querySize = queryCnt.size();
    queryCnt.push_back(INT_MAX);
}

int SimSearcher::jaccardT(double threshold) {
    return max(querySize * threshold,
           (querySize + minSubStrSize) * threshold / (1 + threshold));

}

int SimSearcher::searchJaccard(const char *query, double threshold, vector<pair<unsigned, double> > &result) {
    result.clear();
    newIdx.clear();
    qLen = strlen(query);
    getListsJac(query);
    mergeskip(jaccardT(threshold), INT_MAX);
    for (int i = newIdx.size()-1; i >= 0; --i) {
        double tmpD = calDistJac(newIdx[i], threshold);
        if (tmpD > threshold - EPS)
            result.push_back(make_pair(newIdx[i], tmpD));
    }
    sort(result.begin(), result.end());
    return SUCCESS;
}

int SimSearcher::searchED(const char *query, unsigned threshold, vector<pair<unsigned, unsigned> > &result) {
    result.clear();
    newIdx.clear();
    qLen = strlen(query);
    getListsED(query);
    mergeskip(querySize-threshold*q, threshold);
    int size = smallStr.size();
    for (int i = newIdx.size()-1; i >= 0; --i) {
        int idx = newIdx[i];
        unsigned tmpU = calDistED(dataStr[idx], query, threshold);
        if (tmpU <= threshold)
            result.push_back(make_pair(idx, tmpU));
    }
    for (int j = 0; j < size; ++j) {
        int tmp1 = smallStr[j];
        //if (abs(len[idx]-squerysize)<=threshold)
        {
            unsigned tmp2 = calDistED(dataStr[tmp1], query, threshold);
            if (tmp2 <= threshold)
                result.push_back(make_pair(tmp1, tmp2));
        }
    }
    sort(result.begin(), result.end());
    return SUCCESS;
}

