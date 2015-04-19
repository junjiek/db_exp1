#include "SimSearcher.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstring>
#include <unordered_set>
#include <climits>

using namespace std;

int visitWord[MAXN];
int wordExis[MAXN][BUFSIZE];
int globalWordIdx[MAXN];
int n_Hashq[BUFSIZE];

SimSearcher::SimSearcher() {
    minSubStrSize = INT_MAX;
    maxListSize = 0;
    letterNum = 0;
    wordNum = 0;
    visitor = 0;
    memset(globalWordIdx, -1, sizeof(globalWordIdx));
    dataStr.clear();
    smallStr.clear();
    invertedListED.clear();
}

SimSearcher::~SimSearcher() {}

inline int mypow(int x, int y) {
    int result = 1;
    for (int i = 0; i < y; i ++)
        result *= x;
    return result;
}

void SimSearcher::prepareHash() {
    int Hashq = mypow(HASH, q);
    // n_Hashq[n] = n*(HASH^q)
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
        buildED(buf, lineNum);
        buildJac(buf, lineNum);
    };
    fin.close();
    visitLine.resize(dataStr.size());
    return SUCCESS;

}

void SimSearcher::buildED(const char * str, int lineNum) {
    if (lineLen[lineNum] < q) {
        smallStr.push_back(lineNum);
        return;
    }
    int hashCode = 0;
    for (int i = 0; i < q; i++) {
        hashCode = hashCode * HASH + str[i];
    }
    invertedListED[hashCode].push_back(lineNum);
    for (int i = q; i < lineLen[lineNum]; i++) {
        hashCode = hashCode * HASH - n_Hashq[(int)(str[i-q])] + str[i];
        vector<int> &list = invertedListED[hashCode];
        // gram first appear in str
        if (list.empty() || list.back() != lineNum) {
            list.push_back(lineNum);
        }
    }
}

void SimSearcher::buildJac(const char * str, int lineNum) {
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
            if (invertedListJac[idx].empty() ||
                invertedListJac[idx].back() != lineNum) {
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


double SimSearcher::calDistJac(int index, double thershold) {
    vector<int> &wordIdx = wordIdxJac[index];
    int dataSize = wordIdx.size();
    if (dataSize * thershold > querySize
        || querySize * thershold > dataSize)
        return 0;
    int intersec = 0;
    int i = 0, j = 0;
    // cal intersec num in two sorted vector
    while (i < dataSize) {
        while (wordIdx[i] > queryIdxJac[j]) {
            ++j;
        }
        if (wordIdx[i++] == queryIdxJac[j]) {
            ++intersec;
            ++j;
        }
    }
    return (double)intersec / (dataSize + querySize - intersec);
}

void SimSearcher::getListsED(const char* query) {
    possibleLists.clear();
    querySize = qLen + 1 - q;
    if (qLen < q) return;
    int hashCode = 0;
    for (int i = 0; i < q; i++) {
        hashCode = hashCode * HASH + query[i];
    }
    unordered_map<int, vector<int>>::iterator iter;
    iter = invertedListED.find(hashCode);
    if (iter != invertedListED.end()) {
        possibleLists.push_back(&(iter->second));
        if ((iter->second).size() > maxListSize)
            maxListSize = (iter->second).size();
    }
    for (int i = q; i < qLen; i++) {
        hashCode = hashCode * HASH - n_Hashq[(int)(query[i-q])] + query[i];
        iter = invertedListED.find(hashCode);
        if (iter != invertedListED.end()) {
            possibleLists.push_back(&(iter->second));
            if ((iter->second).size() > maxListSize)
                maxListSize = (iter->second).size();
        }
    }
}
void SimSearcher::getListsJac(const char* query) {
    possibleLists.clear();
    queryIdxJac.clear();
    int newWord = 0;
    int curr = 0;

    visitor++;
    bool find = false;
    for (int i = 0; i <= qLen; ++i) {
        if (i == qLen || query[i] == ' ') {
            int idx = globalWordIdx[curr];
            // the word has appeared in dataset && has not appeared in query
            if (idx != -1 && !find) {
                possibleLists.push_back(&invertedListJac[idx]);
                if (invertedListJac[idx].size() > maxListSize)
                    maxListSize = invertedListJac[idx].size();
                if (visitWord[idx] != visitor) {
                    visitWord[idx] = visitor;
                    queryIdxJac.push_back(idx);
                }
            } else {
                newWord++;
            }
            curr = 0;
            find = false;
        } else {
            if (find)
                continue;
            int next = wordExis[curr][(int)query[i]];
            if (next == 0) {
                find = true;
                continue;
            }
            curr = next;
        }
    }
    sort(queryIdxJac.begin(), queryIdxJac.end());
    querySize = queryIdxJac.size() + newWord;
    queryIdxJac.push_back(INT_MAX);
}

void SimSearcher::sortListLen(int b, int e, int len) {
    int i = b, j = e;
    unsigned pivot = possibleLists[(b + e)/2]->size();
    do {
        while (possibleLists[i]->size() > pivot)
            i++;
        while (possibleLists[j]->size() < pivot)
            j--;
        if (i <= j) {
            if (possibleLists[i]->size() != possibleLists[j]->size())
                swap(possibleLists[i], possibleLists[j]);
            i++;
            j--;
        }
    } while (i <= j);
    if (i <= len)
        sortListLen(i, e, len);
    if (j > len)
        sortListLen(b, j, len);
}

void SimSearcher::divideSkip(int T, int threshold) {
    if (T  == 0) {
        for (int i = 0; i < (int)dataStr.size(); i++) {
            if(abs(lineLen[i] - qLen) <= threshold)
                rawResult.push_back(i);
        }
        return;
    }
    ++visitor;
    int longListNum = T - 1;
    int posListSize = possibleLists.size();
    if (longListNum < posListSize && longListNum > 0)
        sortListLen(0, posListSize - 1, longListNum - 1);
    for (int i = longListNum; i < posListSize; i++) {
        for (int idx : *(possibleLists[i])) {
            if (visitLine[idx] != visitor) {
                visitLine[idx] = visitor;
                if (abs(lineLen[idx] - qLen) <= threshold)
                    rawResult.push_back(idx);
            }
        }
    }
}

int SimSearcher::edT(unsigned threshold) {
    return max(0, querySize - (int)(threshold * q));
}

int SimSearcher::searchED(const char *query, unsigned threshold, vector<pair<unsigned, unsigned> > &result) {
    result.clear();
    rawResult.clear();
    qLen = strlen(query);
    getListsED(query);
    divideSkip(edT(threshold), threshold);
    for (int i = rawResult.size() - 1; i >= 0; i--) {
        int idx = rawResult[i];
        unsigned ed = calDistED(dataStr[idx], query, threshold);
        if (ed <= threshold)
            result.push_back(make_pair(idx, ed));
    }
    for (int j = 0; j < (int)smallStr.size(); j++) {
        int idx = smallStr[j];
        unsigned ed = calDistED(dataStr[idx], query, threshold);
        if (ed <= threshold)
            result.push_back(make_pair(idx, ed));
    }
    sort(result.begin(), result.end());
    return SUCCESS;
}

int SimSearcher::jaccardT(double threshold) {
    return max(querySize * threshold,
               (querySize + minSubStrSize) * threshold / (1 + threshold));

}

int SimSearcher::searchJaccard(const char *query, double threshold, vector<pair<unsigned, double> > &result) {
    result.clear();
    rawResult.clear();
    qLen = strlen(query);
    getListsJac(query);
    divideSkip(jaccardT(threshold), INT_MAX);
    for (int i = rawResult.size() - 1; i >= 0; i--) {
        double jac = calDistJac(rawResult[i], threshold);
        if (jac > threshold - EPS)
            result.push_back(make_pair(rawResult[i], jac));
    }
    sort(result.begin(), result.end());
    return SUCCESS;
}




