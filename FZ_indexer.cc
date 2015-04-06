#include "FZ_indexer.h"

/*
 * This method creates an index on the strings stored in the given data file name.
 * The format of the file is as follows:
	* each line represents a string.
	* The ID of each string is its line number, starting from 0.
 * The index is created and  serialized to a file on the disk, and the file name is indexFilename.
 */
bool FZ_Indexer::CreateIndex(const char * dataFilename, unsigned q, const char * indexFilename) {
	global_dataFilename = dataFilename;
	global_indexFilename = indexFilename;
	global_q = q;
	ifstream fin(dataFilename);
	if (!fin)
	{
		return FAILURE;
	}
	string s,word;
	int index =0;
	while (getline(fin,s))
	{
		count_map res;
		if (s[s.length()-1]=='\r')
			s = s.substr(0,s.length()-1);
		give_str[index] = s;
		if ((index == 0) || (s.length()<smin ))
			smin = s.length();

		int pos = 0;
		while(pos+q <= s.length())
		{
			word = s.substr(pos,q);
			count_map::iterator iter;
			iter = res.find(word);
			if (iter != res.end())
			{
				int t = iter->second;
				iter->second = t + 1;
			}
			else
			{
				res.insert(pair<string,int>(word,1));
			}
			pos++;
		}

		count_map::iterator iter;
		for (iter = res.begin();iter != res.end();++iter)
		{
			hash_map::iterator global_iter;
			global_iter = myindex.find(iter->first);
			if (global_iter != myindex.end())
			{
				(global_iter->second).push_back(pair<int,int>(index,iter->second));
			}
			else
			{
				vector<pair<int,int> > myvector;
				myvector.push_back(pair<int,int>(index,iter->second));
				myindex.insert(pair<string, vector<pair<int,int> > >(iter->first,myvector));
			}
		}
		index++;
	}
    lines = index;

	ofstream fout(indexFilename);
	if (!fout)
	{
		return FAILURE;
	}
	hash_map::iterator global_iter;
	for (global_iter = myindex.begin();global_iter != myindex.end();++global_iter)
	{
		fout << (global_iter->first) << " " << (global_iter->second).size();
		vector<pair<int,int> >::iterator it;
		for (it = (global_iter->second).begin();it != (global_iter->second).end();++it)
		{
			fout <<" " << it->first << " " << it->second;
		}
		fout << endl;
	}
	fin.close();
	fout.close();
    return SUCCESS;
}


/*
 * This method should destroy the index and delete the correspond index file on
 * disk (if any).
 */

bool FZ_Indexer::DestroyIndex() {
	myindex.clear();
	lines = 0;
	char cmd[200];
	sprintf(cmd,"./%s",global_indexFilename);
	int res = std::remove(cmd);
	if (res == 0)
		return SUCCESS;
	else
		return FAILURE;
}

/*
 * This method should load the index from the disk into memory. If it's not in memory.
 * Return an error if the index has not been constructed.
 */

bool FZ_Indexer::LoadIndex() {

	if (!(myindex.empty()))
		return SUCCESS;
	//myindex.clear();
	fstream _file;
	_file.open(global_indexFilename,ios::in);
	if(!_file)
		return FAILURE;
	ifstream fin;
	fin.open(global_indexFilename);
	if (!fin)
	{
		return FAILURE;
	}
	string word;
	int num,index,count;
	while (!fin.eof())
	{
		fin >> word;
		fin >> num;
		vector<pair<int,int> > myvector;
		for (int i = 0; i < num; i++)
		{
			fin >> index >> count;
			myvector.push_back(pair<int,int>(index,count));
		}
		myindex.insert(pair<string, vector<pair<int,int> > >(word,myvector));
	}
	fin.close();
	ifstream fin_data(global_dataFilename);
	if (!fin_data)
	{
		return FAILURE;
	}
	index =0;
	string s;
	while (getline(fin_data,s))
	{
		if (s[s.length()-1]=='\r')
			s = s.substr(0,s.length()-1);
		give_str[index] = s;
		if ((index == 0) || (give_str[index].length()<smin ))
			smin = give_str[index].length();
		index++;
	}
	lines = index;
	fin_data.close();
	return SUCCESS;
}

/*
 * It should do a search using the index by finding all the strings in the data
 * file whose edit distance to the query string is within the threshold. The
 * format of result is a pair of integers which respectively stand for the
 * qualified string ID and the edit distance between the qualified string and
 * query string. All results are stored in a vector, sorted based on the
 * qualified string IDs in an ascending order. Return an error if the index is
 * not constructed or not loaded in memory.
 */

bool FZ_Indexer::SearchED(const char *query, unsigned threshold, vector< pair<unsigned, unsigned> > &results) {
	string s(query);
	string word;
	int MIN_COUNT = s.length()-global_q+1-threshold*global_q;
	int pos = 0;
	count_map res;
	intint_map searchED;
	while(pos+global_q <= s.length())
	{
		word = s.substr(pos,global_q);
		count_map::iterator iter;
		iter = res.find(word);
		if (iter != res.end())
		{
			int t = iter->second;
			iter->second = t + 1;
		}
		else
		{
			res.insert(pair<string,int>(word,1));
		}
		pos++;
	}
	count_map::iterator iter;
	for (iter = res.begin(); iter != res.end();++iter)
	{
		word = iter->first;
		int my_count = iter->second;
		hash_map::iterator global_iter;
		global_iter = myindex.find(word);
		if (global_iter != myindex.end())
		{
			vector<pair<int,int> >::iterator it;
			for (it = (global_iter->second).begin();it != (global_iter->second).end();++it)
			{
				intint_map::iterator searchED_iter;
				searchED_iter = searchED.find(it->first);
				if (searchED_iter != searchED.end())
				{
					if (my_count < it->second)
						searchED_iter->second += my_count;
					else
						searchED_iter->second += it->second;
				}
				else
				{
					if (my_count < it->second)
						searchED.insert(pair<int,int>(it->first,my_count));
					else
						searchED.insert(pair<int,int>(it->first,it->second));
				}
			}
		}
	}
	intint_map::iterator searchED_iter = searchED.begin();
	searchED_iter ++;
	for(;searchED_iter != searchED.end();)
	{
		if (searchED_iter->second < MIN_COUNT)
			searchED.erase(searchED_iter++);
		else
			searchED_iter++;
	}
	intint_map ordered(searchED.begin(), searchED.end());
	unsigned ed[120][120];

	for(searchED_iter = searchED.begin();searchED_iter != searchED.end();++searchED_iter)
	{
		int index = searchED_iter->first;
		string ordi = give_str[index];
		for (int i = 0; i < 120; i++)
		for (int j = 0; j < 120; j++)
			ed[i][j] = i+j;
		for (unsigned i = 0;i < s.length();i++)
		for (unsigned j = 0;j < ordi.length();j++)
		{
			if (i>0) ed[i][j] = ed[i-1][j] + 1;
			if ((j>0) && ((i==0)||(ed[i][j] > (ed[i][j-1]+1)))) ed[i][j] = ed[i][j-1] + 1;
			if ((i>0)&&(j>0))
			{
				if ((s[i]==ordi[j]) && (ed[i][j]>ed[i-1][j-1]))
					ed[i][j] = ed[i-1][j-1];
				if ((s[i]!=ordi[j]) && (ed[i][j]>(ed[i-1][j-1]+1)))
					ed[i][j] = ed[i-1][j-1] + 1;
			}
		}
		if (ed[s.length()-1][ordi.length()-1]<=threshold)
		{
			results.push_back(pair<int,int>(index,ed[s.length()-1][ordi.length()-1]));
		}
	}
	if (results.size()>0)
		return SUCCESS;
	else
		return FAILURE;
}


/*
 * It should do a search using the index by finding all the strings in the data
 * file whose jaccard similarity to the query string is not smaller than the threshold. The
 * format of result is a pair of number which respectively stand for the
 * qualified string ID and the jaccard similarity between the qualified string and
 * query string. All results are stored in a vector, sorted based on the
 * qualified string IDs in an ascending order. Return an error if the index is
 * not constructed or not loaded in memory.
 */

bool FZ_Indexer::SearchJaccard(const char *query, double threshold, vector< pair<unsigned, double> > &results) {
	string s(query);
	string word;
	double MIN_COUNT = (s.length()-global_q+1)*threshold;
	double c = (s.length()-2*global_q+2+smin)*threshold/(1+threshold);
	if (c > MIN_COUNT) MIN_COUNT = c;
	int pos = 0;
	count_map res;
	intint_map searchED;
	while(pos+global_q <= s.length())
	{
		word = s.substr(pos,global_q);
		count_map::iterator iter;
		iter = res.find(word);
		if (iter != res.end())
		{
			int t = iter->second;
			iter->second = t + 1;
		}
		else
		{
			res.insert(pair<string,int>(word,1));
		}
		pos++;
	}
	count_map::iterator iter;
	for (iter = res.begin(); iter != res.end();++iter)
	{
		word = iter->first;
		int my_count = iter->second;
		hash_map::iterator global_iter;
		global_iter = myindex.find(word);
		if (global_iter != myindex.end())
		{
			vector<pair<int,int> >::iterator it;
			for (it = (global_iter->second).begin();it != (global_iter->second).end();++it)
			{
				intint_map::iterator searchED_iter;
				searchED_iter = searchED.find(it->first);
				if (searchED_iter != searchED.end())
				{
					if (my_count < it->second)
						searchED_iter->second += my_count;
					else
						searchED_iter->second += it->second;
				}
				else
				{
					if (my_count < it->second)
						searchED.insert(pair<int,int>(it->first,my_count));
					else
						searchED.insert(pair<int,int>(it->first,it->second));
				}
			}
		}
	}
	intint_map::iterator searchED_iter = searchED.begin();
	searchED_iter ++;
	for(;searchED_iter != searchED.end();)
	{
		if (searchED_iter->second < MIN_COUNT)
			searchED.erase(searchED_iter++);
		else
			searchED_iter++;
	}
	intint_map ordered(searchED.begin(), searchED.end());
	for(searchED_iter = searchED.begin();searchED_iter != searchED.end();++searchED_iter)
	{
		int index = searchED_iter->first;
		int similar = searchED_iter->second;

		string ordi = give_str[index];

        double jacc;
        jacc = ((double)similar)/(ordi.length()+s.length()-2*global_q+2-similar);
        if (jacc >= threshold)
        {
			results.push_back(pair<int,double>(index,jacc));
        }
	}
	if (results.size()>0)
		return SUCCESS;
	else
		return FAILURE;
}
