#include "SimSearcher.h"

#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <climits>

using namespace std;

#define MIN(a,b,c) (a<b?(a<c?a:c):(b<c?b:c))
#define T1 1
#define MAXN 1000000
#define BUFFSIZ 300
#define hashPara 97362

int wordExis[MAXN][BUFFSIZ];
int wordIndx[MAXN];
int visit[MAXN];

typedef unordered_map<int, vector<int> > hash_map;

hash_map hashED;

SimSearcher::SimSearcher() {
    minSubStrSize = INT_MAX;
    wordNum = 0;
    memset(wordIndx, -1, sizeof(wordIndx));
    itemNum = 1;
    dataStr.clear();
    hashED.clear();
    miniStr.clear();
}

SimSearcher::~SimSearcher() {}
void SimSearcher::defsort(int from, int to, int length) {
    int i = from, j = to;
    unsigned k = data[(from + to) >> 1]->size();
    do {
        while (data[i]->size()>k) ++i;
        while (data[j]->size()<k) --j;
        if (i <= j)
        {
            if (data[i]->size() != data[j]->size())
            {
                vector<int> *tmp = data[i];
                data[i] = data[j];
                data[j] = tmp;
            }
            ++i;
            --j;
        }
    } while (i <= j);
    if (j > length) defsort(from, j, length);
    if (length >= i) defsort(i, to, length);
}

void SimSearcher::mergeskip(int T, int thershold, int qSiz) {
    if (T < 1) {
        int j = dataStr.size() - 1;
        for (int i = miniStr.size()-1; i >= 0; --i) {
            for (int k = j; k > miniStr[i]; --k)
                if (abs(inputLen[k]-qSiz) <= thershold) new_index.push_back(k);
            j = miniStr[i] - 1;
        }
        for (int k = j; k >= 0; --k)
            if (abs(inputLen[k]-qSiz)<=thershold) new_index.push_back(k);
        return;
    }
    ++times;
    int occur = T1;
    int len = data.size();
    leave = T - occur;
    if (leave < len && leave > 0) defsort(0, len-1, leave - 1);
    int i = leave;
    while(i < len) {
        vector<int> &current = *(data[i]);
        for (int j = current.size() - 1; j >= 0; --j) {
            int temp = current[j];
            if (visitor[temp] != times) {
                visitor[temp] = times;
                if (abs(inputLen[temp] - qSiz) <= thershold) new_index.push_back(temp);
            }
        }
        i++;
    }
}
double SimSearcher::calJCD(int index, double thershold) {
    vector<int> &temp = indexJCD[index];
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

unsigned SimSearcher::calED(const char *a, int thershold, int asize,int qSiz, const char* query) {
    int length = qSiz;
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

void SimSearcher::createED(int IDLine, const char * s) {
    int length = inputLen[IDLine];
    if (length < q){miniStr.push_back(IDLine);return;}
    int cur1 = 0,i = 0;
    while(i < q) cur1 = cur1 * hashPara + s[i++];
    hashED[cur1].push_back(IDLine);
    for (int i = q; i < length; ++i) {
        cur1 = cur1 * hashPara + s[i] - hash_v[(unsigned)(s[i-q])];
        vector<int> &arr = hashED[cur1];
        if (arr.empty() || arr.back() != IDLine) arr.push_back(IDLine);
    }
}

void SimSearcher::createJCD(int IDLine, const char * s) {
    int ssize = inputLen[IDLine], sum = 0;
    vector<int> iList;
    vector<int> empty;
    iList.clear();
    empty.clear();
    int current = 1,i = 0;
    while(i < ssize) {
        if (s[i] == ' ') {
            int &tmpI = wordIndx[current];
            if (tmpI == -1) {
                tmpI = wordNum++;
                listJCD.push_back(empty);
            }
            vector<int> &arr = listJCD[tmpI];
            if (arr.empty() || arr.back() != IDLine) {
                arr.push_back(IDLine);
                iList.push_back(tmpI);
                ++sum;
            }
            current = 1;
        } else {
            int to = s[i];
            int &next = wordExis[current][to];
            if (next == 0) next = ++itemNum;
            current = next;
        }
        ++i;
    }
    int &tmpI = wordIndx[current];
    if (tmpI == -1) {
        tmpI = wordNum++;
        listJCD.push_back(empty);
    }
    vector<int> &arr = listJCD[tmpI];
    if (arr.empty() || arr.back() != IDLine) {
        arr.push_back(IDLine);
        iList.push_back(tmpI);
        ++sum;
    }
    if (minSubStrSize > sum) minSubStrSize = sum;
    sort(iList.begin(), iList.end());
    indexJCD.push_back(iList);
}

void SimSearcher::hash_init() {
    hash_v[1] = 1;
    for (int i = 0; i < q; ++i) hash_v[1] = hash_v[1] * hashPara;
    hash_v[0] = 0;
    for (int i = 2; i < BUFFSIZ; ++i) hash_v[i] = hash_v[i-1]+hash_v[1];
}

int SimSearcher::createIndex(const char *filename, unsigned q) {
    this-> q = q;
    hash_init();
    ifstream fin(filename);
    string line;
    char * buf;
    while (getline(fin, line)) {
        inputLen.push_back((int)line.length());
        int lineNum = dataStr.size();
        buf = (char*)malloc(1000);
        strcpy(buf, line.c_str()); 
        dataStr.push_back(buf);
        createED(lineNum, buf);
        createJCD(lineNum, buf);
    };
    fin.close();
    visitor.resize(dataStr.size());
    return SUCCESS;

}

void SimSearcher::EDSets(int qSiz, const char* query){
    data.clear();
    querySize = qSiz + 1 - q;
    if (qSiz < q) return;
    int curent = 0;
    for (int i = 0; i < q; ++i) curent = curent * hashPara + query[i];
    hash_map::iterator it = hashED.find(curent);
    if (it != hashED.end()) data.push_back(&(it->second));
    for (int i = q; i < qSiz; ++i) {
        curent = curent * hashPara + query[i] - hash_v[(unsigned)(query[i-q])];
        hash_map::iterator it = hashED.find(curent);
        if (it == hashED.end()) continue;
            data.push_back(&(it->second));
    }
}
void SimSearcher::JCDSets(int qSiz, const char* query) {
    data.clear();
    querySize = 0;
    otherWord = 0;
    queryCnt.clear();
    int current = 1;
    ++times;
    bool find = false;
    //cout << squery << endl;
    for (int i = 0; i < qSiz; ++i)
        if (query[i] == ' ') {
            int num = wordIndx[current];
            if (!find && num != -1)
            {
                data.push_back(&listJCD[num]);
                if (visit[num] != times)
                {
                    visit[num] = times;
                    queryCnt.push_back(num);
                }
                //cout << "push" << tmpI << endl;
            }
            else
                ++otherWord;
            current = 1;
            find = false;
        }
        else
        {
            if (find)
                continue;
            int to = query[i];
            int &next = wordExis[current][to];
            if (next == 0)
            {
                find = true;
                continue;
            }
            current = next;
        }
    {
        int num = wordIndx[current];
        if (!find && num != -1) {
            data.push_back(&listJCD[num]);
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
    int qSiz = strlen(query);
    JCDSets(qSiz,query);

    mergeskip(ceil(fmax(threshold*querySize, (querySize+minSubStrSize)*threshold/(1+threshold))), INT_MAX,qSiz);
    for (int i = new_index.size()-1; i >= 0; --i) {
        double tmpD = calJCD(new_index[i], threshold);
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
    int qSiz = strlen(query);
    // char* query = (char*)query;
    EDSets(qSiz,query);
    mergeskip(querySize-threshold*q, threshold,qSiz);//MERGE
    int size = miniStr.size();
    for (int i = new_index.size()-1; i >= 0; --i) {
        int tmpI = new_index[i];
        unsigned tmpU = calED(dataStr[tmpI], threshold, inputLen[tmpI],qSiz,query);
        if (tmpU <= threshold)
            result.push_back(make_pair(tmpI, tmpU));
    }
    for (int j = 0; j < size; ++j) {
        int tmp1 = miniStr[j];
        //if (abs(len[tmpI]-squerysize)<=threshold)
        {
            unsigned tmp2 = calED(dataStr[tmp1], threshold, inputLen[tmp1],qSiz,query);
            if (tmp2 <= threshold)
                result.push_back(make_pair(tmp1, tmp2));
        }
    }
    sort(result.begin(), result.end());
    return SUCCESS;
}

