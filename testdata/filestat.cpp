#include <fstream>
#include <limits.h>
#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;

#define MAX_ALPH_SIZE (UCHAR_MAX+1)

void printUsage( const char *cmd ) {
	cerr << "Usage: " << cmd << " FILE" << endl;
	cerr << "Requires a file and outputs its file name, alphabet size, number of newlines and size" << endl;
}

int main(int argc, char **argv) {
	if (argc != 2) {
		printUsage( argv[0] );
		return 1;
	}
	ifstream in( argv[1] );
	if (in.fail()) {
		printUsage( argv[0] );
		cerr << endl << "Unable to open file " << argv[1] << endl;
		return 1;
	}

	//collect alphabet information and file size
	long alph[MAX_ALPH_SIZE];
	for (int i = 0; i < MAX_ALPH_SIZE; i++) alph[i] = 0;
	int alphsize = 0;
	long n = 0;
	for (int c = in.get(); c != EOF; c = in.get()) {
		if (alph[c] == 0) {
			alphsize++;
		}
		alph[c]++;
		n++;
	}

	//and print results
	cout << alphsize << " " << (alph[(int)'\n']+1) << " " << n;
	return 0;
}
