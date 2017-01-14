//////////////////////////////////////////////////////////////////////////////*
// Herman Stehouwer
// j.h.stehouwer@uvt.nl
////////////////////////////////////////////////////////////////////////////////
// Filename: main.cpp
// Licenced under the GPLv3, see the LICENCE file.
//
// Copyright (C) 2010 Herman Stehouwer
// //
// // This program is free software: you can redistribute it and/or modify
// // it under the terms of the GNU General Public License as published by
// // the Free Software Foundation, either version 3 of the License, or
// // (at your option) any later version.
// //
// // This program is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// // GNU General Public License for more details.
// //
// // You should have received a copy of the GNU General Public License
// // along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
// This file is part of the suffix array package.
////////////////////////////////////////////////////////////////////////////////
// This file contains two simple examples of using the suffixarray package.
// One with a simple string (not commented)
// And one with a list of integers read in from a file.
//////////////////////////////////////////////////////////////////////////////*/

#include <sstream>
#include <iostream>
#include <string>
#include "intstring.h"
#include "wordstring.h"
#include <vector>
#include "suffixarray.h"
#include <fstream>
#include <cstdlib>
#include <limits>
#include <math.h>
#include <map>

#define HAVE_GETOPT_H

#ifdef HAVE_GETOPT_H
#define GNU_SOURCE

#include <getopt.h>

#else
extern "C" {
  string optarg;
  extern int optind, opterr, optopt;
  struct option { const string name; int has_arg; int *flag; int val; };
#define no_argument            0
#define required_argument      1
#define optional_argument      2
#ifdef HAVE_GETOPT_LONG_ONLY
  extern int getopt_long_only (int argc, string const argv[],
    const string optstring, const struct option *longopts, int
    *longindex);
#else
#warning \
Gnu Getopt Library not found: \
cannot implement long option handling

  extern int getopt(int argc, string const argv[], const string
      optstring);
  inline int getopt_long_only(int argc, string const argv[], const char
      *optstring, const struct option *longopts, int *longindex) {
    return getopt(argc, argv, optstring);
  }
#endif
} // extern "C"
#endif // #ifdef HAVE_GETOPT_H #else

static struct option long_options[] = {
        {"help", no_argument,       0, 'h'},
        {"file", required_argument, 0, 'f'},
        {0,      0,                 0, 0}
};

using namespace std;
using namespace ns_suffixarray;

typedef suffixarray <wordstring> sa_string;
typedef suffixarray <intstring> sa_int;
typedef suffixarray <string> sasa;

string program_name;

void usage() {
    cerr << "Usage: " << program_name << "[OPTION]..." << endl;
    cerr << "This program reads in a corpus and stores it in a ";
    cerr << "suffixarray.  It then returns the number of occurrences ";
    cerr << "of strings in the corpus." << endl;
    cerr << "  -h, --help        ";
    cerr << "Show this help and exit" << endl;
    cerr << "  -f, --file FILE   ";
    cerr << "Filename of the corpus to be read" << endl;
    cerr << "  -s, --server PORT ";
    cerr << "Turn server mode on, listening on the port" << endl;
    exit(0);
}

class Mismatcher {
protected:
    string text;
    string pattern;
    int textLength;
    int patternLength;
    int k;
    int* subsetPositions;
    int subsetCount;

    Mismatcher(string text, int textLength, string pattern, int patternLength, int kVal);

    void defineSubsets(int every);

    // This is abstract method, this needs to be implemented
    virtual int updateMismatches(int positionInText) = 0;
    virtual void init() = 0;

public:
    int* findMismatches();

};

Mismatcher::Mismatcher(string t, int tLength, string p, int pLength, int kVal) {
    text = t;
    textLength = tLength;
    pattern = p;
    patternLength = pLength;
    k = kVal;
}

void Mismatcher::defineSubsets(int every) {
    subsetCount = 0;
    for (int i = 0; i < textLength; i += every) {
        subsetCount++;
    }

    subsetPositions = new int[subsetCount];
    for(int c = 0; c < subsetCount; c++) subsetPositions[c] = 0;

    int counter = 0;
    for (int i = 0; i < textLength; i += every) {
        subsetPositions[counter++] = i;
    }
}

int* Mismatcher::findMismatches() {
    int* mismatches = new int[textLength];
    // non checked mismatches will have -1
    for(int c = 0; c < textLength; c++) mismatches[c] = -1;

    defineSubsets(1);

    int i = 0;

    while (i < textLength) {
        int l = 0;
        int j = 0;
        int targetJ = 0;

        while (j + l < patternLength) {
            int tmpL = 0;
            int tmpI = i;
            int tmpJ = j;
            while (tmpJ < patternLength && text[tmpI++] == pattern[tmpJ++]) {
                tmpL++;
            }

            if (tmpL > l) {
                l = tmpL;
                targetJ = j;
            }
            j++;
        }

        for (int sCounter = 0; sCounter < subsetCount; sCounter++) {
            // position unset is done with -1
            if (subsetPositions[sCounter] < 0) {
                continue;
            }

            int pos = subsetPositions[sCounter];
            if (pos <= i && i < pos + patternLength) {
                int mismatchAmount = updateMismatches(pos);
                if (mismatchAmount > k) {
                    subsetPositions[sCounter] = -1;
                }

                if (mismatchAmount > mismatches[pos]) {
                    mismatches[pos] = mismatchAmount;
                }
            }

        }

        if(l > 0) {
            i += l;
        } else {
            i += 1;
        }
    }

    return mismatches;
}


/**
 * Derived class
 */
class NaiveMismatcher : public Mismatcher {
protected:
    int updateMismatches(int positionInText);
    void init();

public:
    NaiveMismatcher(string text, int textLength, string pattern, int patternLength, int kVal);
};

NaiveMismatcher::NaiveMismatcher(string text, int textLength, string pattern, int patternLength, int kVal) :
        Mismatcher(text, textLength, pattern, patternLength, kVal) {
    init();
};

void NaiveMismatcher::init() {
    // some empty space initialization
    cout << "Naive mismatcher init" << endl;
}

int NaiveMismatcher::updateMismatches(int positionInText) {
    int mismatches = 0;
    for (int i = 0; i < patternLength; i++) {
        if (text[positionInText + i] != pattern[i]) {
            mismatches++;
        }
        if (mismatches > k) {
            return mismatches;
        }
    }

    return mismatches;
}

class KangarooMismatcher : public Mismatcher {

protected:
    vector<string> suffixArray;
    map<string, int> lcpIndices;
    vector<int> lcp;
    vector<int> lcpGroups;
    int groupSize;

    int updateMismatches(int positionInText);

    void init();

    void createLCP();

    int findIndex(string text);

    int RMQ(string text1, string text2);

public:
    KangarooMismatcher(vector<string> suffixArray, string text, int textLength, string pattern, int patternLength, int kVal);
};

KangarooMismatcher::KangarooMismatcher(vector<string> suffixArray, string text, int textLength, string pattern, int patternLength, int kVal)
        : Mismatcher(text, textLength, pattern, patternLength, kVal) {
    this->suffixArray = suffixArray;
    init();
}

void KangarooMismatcher::init() {
    this->groupSize = (int) round(sqrt(suffixArray.size()));
    createLCP();
}

void KangarooMismatcher::init() {
    createLCP();
}

void KangarooMismatcher::createLCP() {
    int lcpGroupMin = 100000;
    for (int i = 0; i < suffixArray.size() - 1; ++i) {

        //Calculate LCP between string and string i + 1
        string string1 = suffixArray[i];
        string string2 = suffixArray[i + 1];
        int string1Length = (int) string1.length();
        int string2Length = (int) string2.length();

        //Compare them for LCP
        int lcpValue = 0;
        int j = 0;
        while (true) {

            if (j == string1Length || j == string2Length) {
                //Write LCP value and set group min
                lcp.push_back(lcpValue);
                lcpIndices.insert(pair<string, int>(string1, i));
                if (lcpValue < lcpGroupMin) {
                    lcpGroupMin = lcpValue;
                }
                break;
            }

            char c1 = string1[j];
            char c2 = string2[j];
            if (c1 == c2) {
                lcpValue++;
            } else {
                //Write LCP value and set group min
                lcp.push_back(lcpValue);
                lcpIndices.insert(pair<string, int>(string1, i));
                if (lcpValue < lcpGroupMin) {
                    lcpGroupMin = lcpValue;
                }
                break;
            }

            j++;
        }

        //Find group minimum
        if ((i + 1) % groupSize == 0) {
            if ((i + 1) / groupSize > 0) {
                lcpGroups.push_back(lcpGroupMin);
            }
            lcpGroupMin = 1000000;
        }
    }
    lcp.push_back(0);
    lcpIndices.insert(pair<string, int>(suffixArray[suffixArray.size() - 1], (const int &) (suffixArray.size() - 1)));

    // Some weird edge case
    if ((suffixArray.size() - 1) % groupSize != 0) {
        if (lcpGroupMin == 1000000) {
            lcpGroupMin = 0;
        }
        lcpGroups.push_back(lcpGroupMin);
    }
}

int KangarooMismatcher::findIndex(string text) {
    return lcpIndices.find(text)->second;
}

int KangarooMismatcher::RMQ(string text1, string text2) {
    int lcpIndex1 = findIndex(text1);
    int lcpIndex2 = findIndex(text2);

    if (lcpIndex1 == lcpIndex2) return (int) text1.length();
    if (lcpIndex1 > lcpIndex2) {
        int tmp = lcpIndex1;
        lcpIndex1 = lcpIndex2;
        lcpIndex2 = tmp;
    }

    int groupIndex1 = lcpIndex1 / groupSize;
    int groupIndex2 = lcpIndex2 / groupSize;

    int rmq = 100000;
    for (int i = groupIndex1 + 1; i <= groupIndex2 - 1; i++) {
        if (lcpGroups[i] < rmq) {
            rmq = lcpGroups[i];
        }
    }
    if (rmq == 0) {
        return rmq;
    }

    if (groupIndex1 == groupIndex2) {
        for (int i = lcpIndex1; i < lcpIndex2; ++i) {
            if (lcp[i] < rmq) {
                rmq = lcp[i];
            }
        }
    } else {
        for (int i = lcpIndex1; i < (groupIndex1 + 1) * groupSize; ++i) {
            if (lcp[i] < rmq) {
                rmq = lcp[i];
            }
        }
        for (int i = lcpIndex2 - 1; i >= groupIndex2 * groupSize; --i) {
            if (lcp[i] < rmq) {
                rmq = lcp[i];
            }
        }
    }
    return rmq;
}


int KangarooMismatcher::updateMismatches(int positionInText) {
    if (positionInText <= textLength - patternLength) {
        int mismatches = 0;
        int l;
        string substring = text.substr(positionInText, patternLength);
        string tempPattern = pattern;

        while (true) {
            l = RMQ(substring, tempPattern) + 1;
            if (l > substring.length()) {
                return mismatches;
            } else {
                mismatches++;
            }
            substring = substring.substr(l, substring.length());
            tempPattern = tempPattern.substr(l, tempPattern.length());
            if (substring.empty() || mismatches > k) {
                return mismatches;
            }
        }
    } else {
        return k + 1;
    }
}

int
main(int argc, string argv[]) {

    char text[] = "Now, you will have a second problem, is that str isn't large enough to hold str2. So you will need to increase the length of it. Otherwise, you will overrun str - which is also undefined behavior.";
    char pattern[] = "overa";
    int tLength = 196;
    int pLength = 5;
    int k = 3;
    //TODO I need suffix array here :)

    NaiveMismatcher mismatcher(text, tLength, pattern, pLength, k);

    int* mismatches = mismatcher.findMismatches();

    for (int i = 0; i < tLength - pLength + 1; i++) {
        if(0 <= mismatches[i] && mismatches[i] <= k) {
            cout << i << ":" << mismatches[i] << endl;
        }
    }

    vector<string> suffixArray = {"a", "an", "ana", "ban", "n", "na", "nan"};

    KangarooMismatcher kangarooMismatcher(suffixArray, "banana", 6, "ana", 3, 1);
    int* kangarooMismatches = kangarooMismatcher.findMismatches();
    for (int i = 0; i < 4; i++) {
        if(0 <= kangarooMismatches[i] && kangarooMismatches[i] <= k) {
            cout << i << ":" << mismatches[i] << endl;
        }
    }

    cout << "Test" << endl;
    return 0;
}


//int
//main(int argc, string argv[]) {
//  program_name = argv[0];
//
//  ifstream is;
//
//  // Handle arguments
//  int opt;
//  int option_index;
//  const string optstring="hf:s:";
//  while ((opt = getopt_long_only(argc, argv, optstring, long_options,
//          &option_index)) !=-1){
//    switch (opt) {
//      case 'h':
//        usage();
//        break;
//      case 'f':
//        is.open (optarg, ios::in );
//        if (!is.good()) {
//          cerr << "cannot open input file " << optarg << endl;
//          exit (-1);
//        }
//        break;
//      default:
//        cerr << "unknown argument " << opt << endl;
//        exit (-1);
//    }
//  }
//  if (!is.is_open()) {
//    cerr << "input file is not open" << endl;
//    exit (-1);
//  }
//
//
//    wordstring word;
//    is >> word;
//    //word.push_back(numeric_limits<int>::max());
//    word.push_back("~~~~~~~~~~~~~");
//    sa_string tree(word);
//    is.close();
//
//  sa_string::size_type result;
//    cerr << "STREE: SUFFIXARRAY Ready!" << endl;
//    char input[1000000];
//    while ( cin.getline(input, 1000000) ) {  // read the string and write it back
//      //cerr << "STree: RECIEVED ***" << input << "***"<< endl;
//    wordstring totest;
//      istringstream a (string(input), istringstream::in);
//      a >> totest;
//      result = tree.find_all_positions_count(totest);
//      //cerr << "SARRAY: Sending result: ***"<<result <<"***"<< endl;
//      cout << result << endl;
//    }
//}
// end of file: main.cpp
