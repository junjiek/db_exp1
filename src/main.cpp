#include "SimSearcher.h"
#include <iostream>
#include <cstdio>

using namespace std;

int main(int argc, char **argv)
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s data query\n", argv[0]);
		return 1;
	}

	SimSearcher searcher;

	vector<pair<unsigned, unsigned> > resultED;
	vector<pair<unsigned, double> > resultJaccard;

	unsigned q = 2, edThreshold = 2;
	double jaccardThreshold = 0.85;

	searcher.createIndex(argv[1], q);
	
	vector<string> strs;
	FILE *fp = fopen(argv[1], "r");
	char *line = NULL;
	size_t len = 0;
	while (getline(&line, &len, fp) != -1) {
		line[strlen(line)-1] = '\0';
		strs.push_back(line);
	}
	fclose(fp);

	fp = fopen(argv[2], "r");
	clock_t start, end;
	while (getline(&line, &len, fp) != -1) {
		line[strlen(line)-1] = '\0';
		// ============= Jaccard =============
		start = clock();
		searcher.searchJaccard(line, jaccardThreshold, resultJaccard);
		end = clock();
		printf("+ %lf %s (%.2lfs passed, %lu results found)\n",
				jaccardThreshold, line, (end-start)/1000.0, resultJaccard.size());
		for (auto s: resultJaccard)
			printf("  %lf %s\n", s.second, strs[s.first].c_str());

		// =========== Edit Distance ===========
		start = clock();
		searcher.searchED(line, edThreshold, resultED);
		end = clock();
		printf("+ %d %s (%.2lfs passed, %lu results found)\n",
			edThreshold, line, (end-start)/1000.0, resultED.size());
		for (auto s: resultED)
			printf("  %d %s\n", s.second, strs[s.first].c_str());
	}
	fclose(fp);

	// ifstream fin(argv[1]);
	// string s;
	// int i = 0;
	// while(getline(fin, s))
	// {
	// 	searcher.searchED(s.c_str(), edThreshold, resultED);
	// 	cerr << s << endl;
	// 	cerr << i ++ << endl;
	// }
	// searcher.searchED("qu", edThreshold, resultED);

	return 0;
}

