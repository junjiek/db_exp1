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

void SimSearcher::defsort(int from, int to, int length) {
    int i = from, j = to;
    unsigned k = rawResult[(from + to) >> 1]->size();
    do {
        while (rawResult[i]->size()>k) ++i;
        while (rawResult[j]->size()<k) --j;
        if (i <= j) {
            if (rawResult[i]->size() != rawResult[j]->size()) {
                vector<int> *tmp = rawResult[i];
                rawResult[i] = rawResult[j];
                rawResult[j] = tmp;
            }
            ++i;
            --j;
        }
    } while (i <= j);
    if (j > length) defsort(from, j, length);
    if (length >= i) defsort(i, to, length);
}

void SimSearcher::mergeskip(int T, int thershold, int qLen) {
    if (T < 1) {
        int j = dataStr.size() - 1;
        for (int i = smallStr.size()-1; i >= 0; --i) {
            for (int k = j; k > smallStr[i]; --k)
                if (abs(lineLen[k]-qLen) <= thershold) new_index.push_back(k);
            j = smallStr[i] - 1;
        }
        for (int k = j; k >= 0; --k)
            if (abs(lineLen[k]-qLen)<=thershold) new_index.push_back(k);
        return;
    }
    ++times;
    int occur = T1;
    int len = rawResult.size();
    leave = T - occur;
    if (leave < len && leave > 0) defsort(0, len-1, leave - 1);
    int i = leave;
    while(i < len) {
        vector<int> &curr = *(rawResult[i]);
        for (int j = curr.size() - 1; j >= 0; --j) {
            int temp = curr[j];
            if (visitor[temp] != times) {
                visitor[temp] = times;
                if (abs(lineLen[temp] - qLen) <= thershold) new_index.push_back(temp);
            }
        }
        i++;
    }
}
double SimSearcher::calJac(int index, double thershold) {
    vector<int> &temp = wordIdxJac[index];
    int length = temp.size(), bsize = querySize;
    if (length * thershold > bsize || bsize * thershold > length) return 0;
    int intersec = 0, q = otherWord + length;
    int i= 0,j = 0;
    while(i < length) {
        while (temp[i] > queryCnt[j]) ++j, ++q;
        if (temp[i++] == queryCnt[j]) ++intersec, ++j;
    }
    return (double)intersec / (q + bsize-j);
}

unsigned SimSearcher::calED(const char *a, int thershold, int asize,int qLen, const char* query) {
    int length = qLen;
    if (abs(asize-length)>thershold) return thershold+1;
    int mm[1024];
    int offset = thershold+1;
    int lo = 1, hi = 2 * thershold + 1;
    int anspos = length - asize + offset;
    int i = 0;
    while(i <= thershold + 1) {mm[offset - i] = mm[offset + i] = i;i++;}
    for (int i = 0; i < asize; ++i) {
        for (int j = lo; j <= hi; ++j) {
            int x = mm[j-1] + 1;
            int y = mm[j];
            int z = mm[j+1] + 1;
            int pos = i + j - offset;
            if (pos < 0 || pos >= length || a[i] != query[pos]) ++y;
            mm[j] = MIN(x,y,z);
        }
        while (mm[lo] > thershold) ++lo;
        while (mm[hi] > thershold) --hi;
        if (lo > hi) return thershold+1;
    }
    return mm[anspos];
}

void SimSearcher::createED(int lineNum, const char * str) {
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

void SimSearcher::createJac(int lineNum, const char * str) {
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
        createED(lineNum, buf);
        createJac(lineNum, buf);
    };
    fin.close();
    visitor.resize(dataStr.size());
    return SUCCESS;

}

void SimSearcher::getListsED(int qLen, const char* query) {
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
void SimSearcher::getListsJac(int qLen, const char* query) {
    rawResult.clear();
    querySize = 0;
    otherWord = 0;
    queryCnt.clear();
    int curr = 0;
    ++times;
    bool find = false;
    //cout << squery << endl;
    for (int i = 0; i < qLen; ++i) {
        if (query[i] == ' ') {
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
    {
        int num = globalWordIdx[curr];
        if (!find && num != -1) {
            rawResult.push_back(&invertedListJac[num]);
            if (visit[num] != times) {
                visit[num] = times;
                queryCnt.push_back(num);
            }
        }
        else
            ++otherWord;
    }
    sort(queryCnt.begin(), queryCnt.end());
    querySize = queryCnt.size();
    queryCnt.push_back(INT_MAX);
}

int SimSearcher::searchJaccard(const char *query, double threshold, vector<pair<unsigned, double> > &result) {
    //cout << "jaccard" << endl;
    result.clear();
    new_index.clear();
    int qLen = strlen(query);
    getListsJac(qLen,query);

    mergeskip(ceil(fmax(threshold*querySize, (querySize+minSubStrSize)*threshold/(1+threshold))), INT_MAX,qLen);
    for (int i = new_index.size()-1; i >= 0; --i) {
        double tmpD = calJac(new_index[i], threshold);
        //cout << filter[i] << ":" << tmpD << endl;
        if (tmpD > threshold - 1e-8)
            result.push_back(make_pair(new_index[i], tmpD));
    }
    sort(result.begin(), result.end());
    return SUCCESS;
}

int SimSearcher::searchED(const char *query, unsigned threshold, vector<pair<unsigned, unsigned> > &result) {
    //cout << "ed" << endl;
    result.clear();
    new_index.clear();
    int qLen = strlen(query);
    // char* query = (char*)query;
    getListsED(qLen, query);
    mergeskip(querySize-threshold*q, threshold,qLen);//MERGE
    int size = smallStr.size();
    for (int i = new_index.size()-1; i >= 0; --i) {
        int idx = new_index[i];
        unsigned tmpU = calED(dataStr[idx], threshold, lineLen[idx],qLen,query);
        if (tmpU <= threshold)
            result.push_back(make_pair(idx, tmpU));
    }
    for (int j = 0; j < size; ++j) {
        int tmp1 = smallStr[j];
        //if (abs(len[idx]-squerysize)<=threshold)
        {
            unsigned tmp2 = calED(dataStr[tmp1], threshold, lineLen[tmp1],qLen,query);
            if (tmp2 <= threshold)
                result.push_back(make_pair(tmp1, tmp2));
        }
    }
    sort(result.begin(), result.end());
    return SUCCESS;
}

