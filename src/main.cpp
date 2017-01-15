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
// One with a simple char* (not commented)
// And one with a list of integers read in from a file.
//////////////////////////////////////////////////////////////////////////////*/

#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <math.h>
#include <map>
#include <set>

#define HAVE_GETOPT_H

#ifdef HAVE_GETOPT_H
#define GNU_SOURCE

#include <getopt.h>
#include <unordered_map>

#else
extern "C" {
  char* optarg;
  extern int optind, opterr, optopt;
  struct option { const char* name; int has_arg; int *flag; int val; };
#define no_argument            0
#define required_argument      1
#define optional_argument      2
#ifdef HAVE_GETOPT_LONG_ONLY
  extern int getopt_long_only (int argc, char* const argv[],
    const char* optchar*, const struct option *longopts, int
    *longindex);
#else
#warning \
Gnu Getopt Library not found: \
cannot implement long option handling

  extern int getopt(int argc, char* const argv[], const char*
      optchar*);
  inline int getopt_long_only(int argc, char* const argv[], const char
      *optchar*, const struct option *longopts, int *longindex) {
    return getopt(argc, argv, optchar*);
  }
#endif
} // extern "C"
#endif // #ifdef HAVE_GETOPT_H #else

static struct option long_options[] = {
        {"help", no_argument,       0, 'h'},
        {"file", required_argument, 0, 'f'},
        {0, 0,                      0, 0}
};

using namespace std;

char *program_name;

void usage() {
    cerr << "Usage: " << program_name << "[OPTION]..." << endl;
    cerr << "This program reads in a corpus and stores it in a ";
    cerr << "suffixarray.  It then returns the number of occurrences ";
    cerr << "of char*s in the corpus." << endl;
    cerr << "  -h, --help        ";
    cerr << "Show this help and exit" << endl;
    cerr << "  -f, --file FILE   ";
    cerr << "Filename of the corpus to be read" << endl;
    cerr << "  -s, --server PORT ";
    cerr << "Turn server mode on, listening on the port" << endl;
    exit(0);
}

class StringPointer {
public:
    char *text;
    int start;
    int end;
    unsigned long hash;

    bool operator<(const StringPointer &rhs) const {
        if (this->hash == rhs.hash) return false;
        int length1 = rhs.length();
        int length2 = this->length();
        if (length1 > length2) {
            for (int i = 0; i < length2; ++i) {
                char c1 = rhs.text[rhs.start + i];
                char c2 = this->text[start + i];
                if (c1 != c2) {
                    return c2 < c1;
                }
            }
            return length2 < length1;
        } else if (length2 > length1) {
            for (int i = 0; i < length1; ++i) {
                char c1 = rhs.text[rhs.start + i];
                char c2 = this->text[start + i];
                if (c1 != c2) {
                    return c2 < c1;
                }
            }
            return length2 < length1;
        } else {
            for (int i = 0; i < length1; ++i) {
                char c1 = rhs.text[rhs.start + i];
                char c2 = this->text[start + i];
                if (c1 != c2) {
                    return c2 < c1;
                }
            }
            return false;
        }
    }

    StringPointer(char *text, int start, int end) {
        this->text = text;
        this->start = start;
        this->end = end;
        this->hash = getHash();
    }

    int length() const {
        return end - start;
    }

    unsigned long getHash() {
        unsigned long h = 5381;
        for (int i = start; i < end; ++i) {
            char c = text[i];
            h= ((h << 5) + h) + c;
        }
        return h;
    }
};

class SuffixArray {
public:
    vector<StringPointer> subsets;

public:
    SuffixArray(char *text, char *pattern, int textLength, int patternLength) {
        set<StringPointer> uniqueSubsets;
        for (int j = 0; j < textLength; j++) {
            uniqueSubsets.insert(StringPointer(text, j, textLength));
        }
        for (int i = 0; i < patternLength; i++) {
            uniqueSubsets.insert(StringPointer(pattern, i, patternLength));
        }

        for (set<StringPointer>::iterator it = uniqueSubsets.begin(); it != uniqueSubsets.end(); it++) {
            subsets.push_back(*it);
        }
    }
};

class Mismatcher {
protected:
    char *text;
    char *pattern;
    int textLength;
    int patternLength;
    int k;
    int *subsetPositions;
    int subsetCount;

    Mismatcher(char *text, int textLength, char *pattern, int patternLength, int kVal);

    void defineSubsets(int every);

    // This is abstract method, this needs to be implemented
    virtual int updateMismatches(int positionInText) = 0;

    virtual void init() = 0;

public:
    int *findMismatches();
};

Mismatcher::Mismatcher(char *t, int tLength, char *p, int pLength, int kVal) {
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
    for (int c = 0; c < subsetCount; c++) subsetPositions[c] = 0;

    int counter = 0;
    for (int i = 0; i < textLength; i += every) {
        subsetPositions[counter++] = i;
    }
}

int *Mismatcher::findMismatches() {
    int *mismatches = new int[textLength];
    // non checked mismatches will have -1
    for (int c = 0; c < textLength; c++) mismatches[c] = -1;

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

        if (l > 0) {
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
    NaiveMismatcher(char *text, int textLength, char *pattern, int patternLength, int kVal);
};

NaiveMismatcher::NaiveMismatcher(char *text, int textLength, char *pattern, int patternLength, int kVal) :
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
//        if (mismatches > k) {
//            return mismatches;
//        }
    }

    return mismatches;
}

class KangarooMismatcher : public Mismatcher {

protected:
    vector<StringPointer> suffixArray;
    unordered_map<unsigned long, long> lcpIndices;
    int* lcp;

    int updateMismatches(int positionInText);

    void init();

    void createLCP();

    long findIndex(StringPointer text);

    int RMQ(StringPointer text1, StringPointer text2);

public:
    KangarooMismatcher(vector<StringPointer> suffixArray, char *text, int textLength, char *pattern, int patternLength, int kVal);
};

KangarooMismatcher::KangarooMismatcher(vector<StringPointer> suffixArray, char *text, int textLength, char *pattern, int patternLength,
                                       int kVal)
        : Mismatcher(text, textLength, pattern, patternLength, kVal) {
    this->suffixArray = suffixArray;
    this->lcp = new int[suffixArray.size()];
    init();
}

void KangarooMismatcher::init() {
    cout << "\nInit Kangaroo mismatcher\n";
    createLCP();
}

void KangarooMismatcher::createLCP() {
    for (int i = 0; i < suffixArray.size() - 1; ++i) {

        //Calculate LCP between string and string i + 1
        StringPointer text1 = suffixArray[i];
        StringPointer text2 = suffixArray[i + 1];
        int string1Length = text1.length();
        int string2Length = text2.length();

        //Compare them for LCP
        int lcpValue = 0;
        int j = 0;
        while (true) {

            if (j == string1Length || j == string2Length) {
                //Write LCP value and set group min
                lcp[i] = lcpValue;
                lcpIndices.insert(pair<unsigned long, int>(text1.hash, i));
                break;
            }

            char c1 = text1.text[text1.start + j];
            char c2 = text2.text[text2.start + j];
            if (c1 == c2) {
                lcpValue++;
            } else {
                //Write LCP value and set group min
                lcp[i] = lcpValue;
                lcpIndices.insert(pair<unsigned long, int>(text1.hash, i));
                break;
            }

            j++;
        }
    }
    lcp[suffixArray.size() - 1] = 0;
    lcpIndices.insert(pair<unsigned long, int>(suffixArray[suffixArray.size() - 1].hash, (const int &) (suffixArray.size() - 1)));
}

long KangarooMismatcher::findIndex(StringPointer text) {
    return lcpIndices[text.getHash()];
}

int KangarooMismatcher::RMQ(StringPointer text1, StringPointer text2) {
    long lcpIndex1 = lcpIndices[text1.getHash()];
    long lcpIndex2 = lcpIndices[text2.getHash()];
    int rmq = 100000;
    if (lcpIndex1 > lcpIndex2) {
        for (int i = lcpIndex2; i < lcpIndex1; ++i) {
            if (lcp[i] < rmq) {
                rmq = lcp[i];
                if (rmq == 0) return 0;
            }
        }
    } else if (lcpIndex2 > lcpIndex1) {
        for (int i = lcpIndex1; i < lcpIndex2; ++i) {
            if (lcp[i] < rmq) {
                rmq = lcp[i];
                if (rmq == 0) return 0;
            }
        }
    }
    return rmq;
}

int KangarooMismatcher::updateMismatches(int positionInText) {
    if (positionInText <= textLength - patternLength) {
        int mismatches = 0;
        int l;
        StringPointer textPointer = StringPointer(text, positionInText, textLength);
        StringPointer patternPointer = StringPointer(pattern, 0, patternLength);

        while (true) {
            l = RMQ(textPointer, patternPointer) + 1;
            patternPointer.start += l;
            if (patternPointer.length() > 0) {
                mismatches++;
                if (mismatches > k) {
                    return mismatches;
                }
            } else if (patternPointer.length() == 0) {
                return mismatches + 1;
            } else {
                return mismatches;
            }
            textPointer.start += l;
        }
    } else {
        return k + 1;
    }
}

int
main(int argc, char *argv[]) {

    char* text = "GACTTATGCCTGATCGCTGTCAGGTCATACGCTTCATTTATGACTTGGCATAACCGGTTTTCTGATGCCACTCGAAGGCACCGTGTTTAACCCTATTCGCGAGGAAATGAGGAACTTCGCCACACCGTGATAACGTTTCTTTATCTTTACCCTGCATCACCCAAGGCTGGCTCGATAATCGAGTTTTGGCGCATCTTCGCATCGGCGTCGTGACAGGGCATCATGCCTACGCTTAGTTAAAATTTGGCTAATGTAGAAATGTTGGCAAGAGAACCGGAAGAGGCGCTTATTACCATCTAAGGATACGTTATAGAAACTTCTCGGTGGATTCATCCTTACGAACAGGGTGCTCGAAGTCGCCCTTCTTATTTCATCTCCTCGAGCACTGGCGATGTGTTTCTCACCTGCTCCGGCTTATTGAACTCCAGCACCGGTTCTGGGGCTGCTCGTTAAAGCCGTTGTTTGCTGGCAAACGGCAAGACCGGTTCCAGGCCTAAAGTTTCCAGTTGGATCACCAGCCTTCAGCGCTGGCGTATCCGCTCGAGCACTTCATTAAGAGGCTTTGCAGTTTCGCGGCATAGTCCGCGGTGCCATCCGTGCCTGGAGCGGCTGTTCGTTCGAGGTCTAATCCCCTTCGCTTTCTTGAAGCGTCTTTGTTGTAATAGAAACGCGTAGGGTGGTCGAGTGGCTGGAGAGTAAGTGGCCCGTTTGCTGTCGGAGTAGGTCAACCTGAAACCGTCGGCACAAACTGCGACTCATCGAACTGAATCCCTTGCCTCTTTAAACACGTCATATACACCGGTTTTAATGGCCTCGACGCCATCATGGTGGCGGTGCCAATTGCCAACCTGCGCAAAATAGCCGGCGCGTGTTGCCGGTACGAAAAATGCGGCAATCCCGGCGCTTAGCAAATTCGCTGTTATGCTGTGGTCTAAATCGTACAATTTGTAATCCGGGTTTTCGGCGTTAAAACGTTGGGCCAGAGAATATCCACCTCTTTACCCAGTTTCCCCTTCCATAGAATTCTCCAGAACCGAATACTTTATGTCACTGCCTGTGCATTCCCCATTGGACGCCAGTCCGATAGCCAGTGCTGAAGCTGTATAATGTAAACGGTTTCATCGTTTTATCTCTCTTGTTGTACCGAATTGCGCGAATTGTCTCCGCGTTTAGCCGCGGGGTAACATGACATGCTCGAACTTACAGAAAAATAACTTTGTTACATTTGTAAGATAGTAAGGTGTCAGAAAGATGACAATAGGCGGTGACGGCGTGGGTGAGGGAAAATGGAGATGGCAACCATGAAAATAGCGAACCATGAATCAAACTCTACATAATTGCTCATCGTTTCATGCCGGATGCGCTAGTACAACGCCTTAAGGCCTGCTATACAAGTACGTGCAAATTCAACATACTTGTCGCCACTCACCCAGTAGGCCTGATAAGCGCAGCGCGCATCAGGCAATTTACATTTGTCATGTCTCAAAAGGAAGTTTTACTCCCTATCAAATCAACGTGTTATTACCCGCTAAATACGCACTTCTCACTAGATTCATTTCGCCATGGATAAGAATAGCATCAGTATCGGAAACCCACTACATTAGGACTTTTCCTCAGCACGATATCGCGATCGCCAGCTTTAGCCGCTCTGCGTTTTGTATGCTCGACGAGACAGAGTCATGCCCTGGTATCGCTCGCAGCTGCTGAGTGGTGTCGAATTGCTGGATGATAATCGGCGCAAGACCGAGCGATGGCTCAGACCCAAGGGCAGTAGCAAACCGCGGGTTAGCGCTCAGGCGCAACGACCAATCGCCAGCATCTGCTGTTCAAGCCGACATGGTGCCGCCCGCTGAATACCGGCGCTCATGCGGCAGACGTGGAAAACAGGCTCATAACCCCACTTATGCGCTCCCTGGAACTGGTCGCGTTCAGCAAAACCCCAAGCCAGGTTCTCTTCCACCGTCATCCGCGAGAAAGGCGCGACAATTCGCTCCACCGCTTCGCGCATGATTTCGCTGTCTGCCAGTCGGTAATGTCTTTATCATCAACAAGAACTTCACGCTGGTGGCACAAATCGTAAGCCAGAAGCTCATGTGCCAGCAAGTGGTTTTCCCCGCGGTTTCGCGCCAATCATGTGGTAAATCGCCCTGATTGATATGCAGGCTCACCTCATTCACAGCGCTCCGATTTGCATGTAGTGGGCGCTGATTTGTCAAAGGACAACATGACTTTTCCATCTTATGCCTCACTAAATAGGGATCGGATTTCACGGGTTATTACGAGGATCACGCTCCGTGTACTCCGTTTGCATAGCGGCGTCCCCTGATTGACCACGTATAAATTCGGTCGGAGAACTTCCCATCACCATGGTTTCATATCGGTGCCAGATCAACAAGATAGTGGTGTTGTATGATTGCGGCATATTTCGGCAATCAGTAGGCCATCGAGCGTTCTCGTGTCTCTTCGGGTTAAGACTTGCCACGCAGGTTCGGTCGAGGCATTAAATCTCCATCTGCGTCACCATGCAGCGGGCACTCACATGAAATGGGGGCCTGTACATTACCAGCGGTACTCGCTGACGGTTGAGGTCCATCCAGCAAACCAATGCGCTCAAGCCAGGTCCCGCGGCCGCGGTCCCGAGCGCTTCGCTGCGGCGCGGGAAAGGAATAACCGTTTGTTCCAAACAGGCCAGAGAAACAGCCCGGTTTTCAGTTGCTGATGCTGCGCCACCAGCAGGTTTGCCAATTACCGTCATTTTTCACGAGTAAAGCACATGCTGAGAAGGTGCGCACCACGCGATGCGGGCAATTTGCTGCCCGGTAAAACCTTCCAGGTGTGCTGCAAGCGCAGTAAAATGGTGCTCATTGAGTTTGTAGAATCCGGTCAGACAGTCAAAATGGCTGTGGTTTTGCCCCTGGAGCCACATGTTGGCCGATCACGGAAAGCATCCCTCCTGCGGGTACAGTTCAAGATTGACGACATGTTGTTCAAGCCAGCAGGCGCCAAGAAGCGCATCAGGCCGTTTAACAGATAATAATGGCTATGACTCATGCTCGCCTGCTCTCCTTTCGCTTGCGCCGTTTGCCATTTGAGCTTGCGGCGCGTCATGGGCAGCAATCCCTGCGACGCCAGAAATTGATGATGATAGGCACCATCAAACCACCGAGCATTAACACTTGGCCTGTATTCGTTGAAATCACGCATCAAGCGCGACACCACCAGCAAACTTGCCGCCAGAATCAAGCGCAAATGCGAGCGATAATAGCCGAGCACCACTATCGCCGCACACAAACCGCCGAATGTTGTAAAACAGGTGAAGGATTCCGGGCTGACAAAGCCCTGAAGCGCCGCAAACAGCGTTTCCGCAAACCGGCAACGAGGCTATTAAAGGCAAGTCAGCTTGATACGACGCGGGCTGCACCCCAGCGAAACGGCAGGCGATTCCGACTTCGCGGCTATTCGCACCCAACCAGCGCGGGAGGTTGAGCATATCAATAGAGGCCTTAACAAGCTCGAGGGACCGAGGGAAAGTACGATCGCGTACTGTTATTAATGCGCAAACAGCAGCCACGGAAACTTGGCTGACATGTTAAATCGAGCGCTTCGAGGCCTGCCATCAATCATGGCAACTGCGGAAATGTCATGTTTACGCGTACTGCTCGTAATTTACATTTCGAGCATCCTGTCCAAATTACATCCGCTGTCACGAGGTGGAGGATAAAATCGCTTACGGAGCAGCACCAAGAAAGACCAATATCCCAGCAGGAGACAGGGTTACATTAATCGAGCGAGTGGCAAGCAGGTCAGGCACCAGCCGTAATGATTACGCAGGTGTGAACCACTAAACGATCGACTCCAATCTCACAGCACCAGCAGGCGCGACTATAGAACAGTTTCGAGCGCCAAGGGAGCTGGAGGATAATGTACAGAAGTTAAACTTCGCCTGGAAAACGACGCTTTGCCACCCTCGTGAAAATCAGCTCGCAAAGCTGACAGACAGAAGGCCTGAGGGATAAATGTTGACTTCAACGCCTTGAAATACATCCGAACCCAGGGTGGCACAATGTATTTCGGTCGCGAGCCCATTTTAATCCATTTCTGGAAAGGCCCATTTAAATCCCGCCCAGCACAAAGGCAAAAGAGCCCTAATGTTACATCACTAGCCACTTAAGCCTGCCACTGGCCCAGCCATGTCGCTCGCGGTTGGTGCCTAATGGGCGACCCAAGTTACCGCTAGGGCACCCACCAGCATGATACCAACCAGGGGCAGCGATTGAGACTTACATCGTTGTCATATTTCTATACCATGCCGATAGGCAAGAAATACCAGGTTTCGCGACTCCGCAACATCCCAGAATCAGGGCGAATCAGCCAATTGACTTCGCTGCTCCGGGGCGAAGGCCCTCCAATCCCCAATGTACAGCCGCGGTAAGGTACCTCATCCCGGAGTATAAATGCCGATGTGTACCCATAAGGGAAACCGGTCCGTATAGAATAATAGCGAATGGCGTACCTGTACTAAGTGCCTGAAGCCATCACAAAGAGTACTCAATCACCCGGGTCCCCGGTCAATGTCTCGTACCGACCGAGCCATTTTCAAATTCTTCGCAGGCGCACAGAGGACTGGAGATGCGGGAATAGCGAATCCCAAATCGTTGCCAGCATCGAGGAAAGGTAACAATCGCAGGTCACGACTTCCGCATGGTAATACAGAGGCAGAGAAGTTTTCGCTATGCCTCCACCATGACCGTCGACCCGAAGGACATTTGGAGCATCTTCAACCTGCCTGAACCTCGTACTATTGTTGCGAGGAAGGTGAGCCGATTGCAGAGATGGAGGTTAGCATCAGGCGCTTGAGTTACGCACCGGGCGGTACATCCACCCGTTCGATACTCCAGCCGGTAGGCGCTGGCAATGACGATTGCGCCGACGAATCCCGCAGCTACCAGCAGCCAGCCGGTATCAATTGCCCATCATCATCAGCGCGGCGATGAGCATAAATGAGACTGTAGCTGCCAATCATATAAAACCTCACGTGGGCGAAGTTGATTAATGCCGATAATGCGAGTGTAAATCGCTCTGTAGCCGATGGCTATCAGCGCGTAATATGCCCAGCGTGACGCTTGTTACAACATCTGCTGCAAGAAATACAAAACTGCTCAGACATAAGGTGATTTCTATAAAACCCGCGAATCTGATTTACGGGCGGTGGAGACCTTAATATTGCTGCCTGTGGATGAACCCTTCGGCGTGCCACTGGAAAGACACCAAATCAAAATCCCTTAAGTACCGCCTTTTCATCCCAGTTTCAGCGGCCCAATCGACATGTTTGCACCGTCAGCTTTAATCTTCACCAGCGCCAGCGGCTCATCGCTGCCGGTAAGCTCAAGGGCATTACCTCAGAGATTGCACCGCTGCGTGAGGTGATCAAGAAAGAGCCGACGGCATCCGATTTCTCTTGTCTGCTTTCAGCGCATCAAAGCTCCTTCGATATTGCCATTACGGGTCATAGCGTTTTGGCATAGTGACCAACATTCCATTTCGGCGATAGCAGCATACCGGCAATGTTCAGCGGTACATTCCGGCCCCATAAAACTGGGGTTTGAGGCCAACGATGAAGCATCTGCCCCATTTCCGGGTAAACTCAGGCCTTAAACGAAGTTAAGTGTTTCTTTTCAGGGCGTAAGCAGTAAGGCGCGGAGAAATCTTTTCTCCCCGGCGGTAATACCGTCGAAGAAGAAAGGACGTTGGCGTTAGCCGTGCCATCAAGTCTTAAGCGAACACGCGCCATACGTGCCATACTGTTGTTTGGTGTAGAATGATGGCCGCGATGCTGGGGGTTCCACCGTCTCAAGAATGTATTTTGCCGCTTGTTGGCCCCTGGGAAGAGTCCAGCCCGGCAAGTACTGGCATAATGTTGATACACGGTGTTGGGTCAGCTCCGGGTTGGTCGCACTCCCAGGCGATCATCAGAATACCTTCGTGTCTTCATAGATATCTGGACGCAGGCTGGGTAGAAGAAGAAAACACAGATGACCAATTAACTGTCATATTAATTCATATTAACGATTGTTGTAGAAGCTCGAAACGGCTTGTTTCGGGTCGCGCATGCGTGACTATTCCAACGTGCCCAAACCAGTTTGACGCCTCTTAAATTCCCCTTTGGCATCTTGATGTCTTTAATTGCCTGAAGGCCTGTTAAATTCCATATCGCCCCACTGGGCAATCGGTGAGCGGACATCGCGCCGACAAACCCTGCGACTTTGTCATGGTCAGCCATAGCGGTGTGAAAATTGCCAGTCACAAATCGTATCCCTGCGATGATAGTTTTCGCATTTCCGTTTCATAGTCAAATCCCATTCGTGATGTTGGTTGTGTTTTATGTTAACAAATCAGACTGTTCTTTATACTGCACTGTTTTGCCTGTCTGATCTTAAGGGGTTAGCGCAGTATTTTGGTCAATAGCGATTAAACCCTATTTTCATAGTCGATTAAGAAACAGATAATATTCTGAAGTCTTACAGAGACTAAACAGAAAATTGCCTTTGTCAGCATAAAATACAACGGCACAAATGAGAAATAATTCACTATCATTCAGGGGATCATGATCTGGACATTTTCATCTCTTCTAATGTTTTAATTTGTAATTATTGCTGTTAAAATTAATCACCTGCCAAAGAAATAAAAGAGAAAGCCTCCGATTAAATTATTTCGCTACACTGGTTTCACTTTGTGATTACACGGGTTACCCATGAAGCTGACCATCATTCGATTGGCGAAAATTGGTAATACAGACCGGGAATAAGGTGCATTAGCAAAGATCACCTGGCCGATGATTCCCCTTCCTCGTTTACAGGTTGACGATAACCACCGTATCTACGCCGCGCGTTTTAACGAGCTTGCGTCGCTAGCCATGCAGGTTGTCCGGGTAACCTTACAAGCGGCACCGAGGGAGCGGGACCTGTAATTCCCTGCGCGTGCGGGAAGTCTACCCGCTTGTCGCGTGTGGGGCATACTCTGCTGAAGAGGTTTGCGTAACAATTCTGAGTTTCATTGCTGGTGGATGGCGGAGTGCAGCGGTGGAAGTCCCCTTGATGAAGGCTTTATGCAGGCGCTGGGGTTTACGGCACAACAGGGCGGCTGGAGAAGTGTCAAATCGTCAAGTTTAACTTCAAAAGTGATATTGCTGATGCGCTACGGGTCTTTATCAGGCGCAATGTGTGTTGCATGTCTACTGATTTCTTTGGATCTGTAGGCCGGATAAGGCGTTTTAATGCCCCACATCCGGCATGAAGCGGACTAGTACTCGATATTAGCAATATTTGCGGCAACCCAAAATTTCGCTTTAATTACCGTAATTTACCTCATCGGTCGCCGTGCTGTTGGCGTGCCAGTCAACCGCTAGAAACTCAATCCAGCTTTTCAGAGGATCGCCTTTTCTCATCCCAGGTCAGCGGTCCCATTACGGTATCCACGCGAGTTCGGTTTTCAGGTATTTGGCGATTTCAGCCGATCGTCATCAAATTTCAGGCCCGCCTGCAAAGATTGCAGCGCGGCGTCAGGTGGTCCAAACGAATCGCGCCACTTGGGTCCTGTTTTACGCTTTGACTGCGTCAACAATCGGGTTTGTTCGCTCGAACCTGCTAACTGCTTAGTTCTTCTATTGGTCACCAGCAGCCCTTCCGCTGATTTGCCCGCAATGTTAGACAGCGAAACGTTGATACACCTTCGCCTCCATAAACTGAGTTTCAGCCCTGCCGCGCGTGCCTGACCGCAGGATTTGCCCCATTTCCGGGTAACCGCCGTGTAGTAAACGAATTACGATATTCTCTTTTCTAAGACCGCGCCACCCGGTGTTGAATCTTTTCCCCGGCCCATCTGCCATCAAAGAAAATACGTTTGCAATTGCCTTTCTTCAGGCCATGTCCTGCAGGCCTCGCGCCAGACCTTCGCCGTATTGCTGTTTGGGCTGAACGTCAGCAATACGTGCAGGTTTCACTTTCTCAAGAATATATTTCGCCGCCGTCAGGCCCCTGGTTACGAGTCCAGGCCGGTGGGTCGCGTACGTACCGGATAGCCACGGGCGGTCAGCTCCGCTGAGGATTGCCGCTGGCGTAAGCTTAAAATCGCCTTCGTGTCTTCGTAGATGTCAGACGCAGGCTTAGTTGATGAAGAACAGAGGTACCAATCACATATTTAAATTGCCATTTATTTAACGACTTTGTTTCACCGCGCAACACGCCTGTTTCGGTCAGGCATCTACAGATAATTTAAGGAACTTGCAGTTTGTTGCCTTTAATGCCGCCTTTAGCGTTGATATCACAACCGCTTATGCGCCGGTAACTCCTGGTCACCGTACTGCGCAACCGGGTGGAAATGGACAATTGCCCCCAAGAACGCGACTTAATATCTTCTGCCAGCCATATTGCTGAATGCCAGCGCGATACATCCTGCCAGTAAACGCGCTTTACCCTTTATGTTCATCCTGGAACCCCATTCTTCTGGTTATTAATTTGTTGTGATGTTGTTGCCATCAATATTTATTTTCGTTTTATGCATGACTACCCGTGCTTTTTAGCAGCATACTTGGCCCAAACATACCGATTTTATGATATTGAAATAGCTATTTGACAGTGTCTATTAACAATCTGCGTGGGGATCAGTTTGCCGGAGGAACTTAATTATTACAGAGGCCCAAAACAAAACCCCGGCCACGCCATCCAGGGTTCTCTGCTTAAGGTAGCGGAAACTTAAGCTTCAATGGCATCAAGACCGCAATTTTCATACCGCGTTTCTTTTGCTCCAGCTGCGTACACGCTCAGCGGAAACGCTTAAAGGTCAGCGGTTCCTCTGCAAACGTGTACTGTTGTCTTCATTTAGCCAGCGCGACTGATGTCATCGGTGGCCGTTGCGTTCGGGATTCCATCTAGTCACATCGCGTCGGTCAGACGGTTTGCTCGCCTGCTCTTCCCAGTTATCATCTTCAATTGGCCGTTGGCAAAGTTAGCTATCGATTGTCATCCTTGGCAGATAGAGCACCGGGAGCCATCGGCTCGTTGAGGAAATCGTGTCGTCGGAAGACAGAGGTCAAAGGTCATGTCTCTGTGCCCGCCATACGTGAATTCCATCAACGATTGTCTTTGCTGGTTACACAGCAGTTCACGGGCCACCATTTCGACTTCATCCTGGTTAAACCAGCCCAGACCGCTGCTTGGTTTACGCAGATTGAAGAAACCAAGTTTGCGCTTGCGCTTTGGTGGTCGCAGGAACTGTTTGTACCGGAGTACTAAGTCCGCAGAACGTATTCGTGGGATCTCTGCTTTGATCACCAGTGAACGGCGAAGGAGACCATCGCACAATTATTCTTCTAGGTTGAAACGGCGCACTGCTTTCATCAGGCCGATGTTACCTTCCTGAATCAAATCGCCTGTGGCAGGCCAGCCCGCATAATTACGAGCATATGAACAACAAACCGCAGGTGAGACAGGCTGCAGCGTTTTAGCTGCTTTTCCGGGCACTGCCCGGTAATGCATATTTTAAGCCATGCCCCGCTCTCCGTCAGCCGACAACATCATACGCGTTGTCGCTGCCCGATGGGAATCCCAGGTTGCCAACTGGCTAAAGCAAACTTGCATATTTGTCAGTCATTCAAATCCTCACGATATCTTCTAGGGCCTGCCTGTCGCAACAAAGCTTGCCAGGATCAAGAGCGAAAGGTTATCATATTCAACTGTTTTATCAGACCAATCTGTTGTATCCACAAGTTCAATTCACCGGATGTGAAATAAATTACGCACAAAATGTGACATAGAGATGAAATACCGGGAGACTGAGGGTCTCTTCCCTGCTACGAACCCAACTTGCAGGGAAAGAGAGTAACACGCTTTATTATTCAGGCTAACCTAGTAAATGTTGTACCGTGGCAGCCACGCTTCGAGCCCAGCCAATCATCGAGCATACCAGCAGCAATAGCAGGGCATTGCCTATCGAATACCGATAAGGCCCATATTGATATCAAACTAACTTCGTTTCCGAGAAAACCTGTGCCACACCTCTTCGCAACCGCCGATGACAATCGCAGCACCAGAATTTTCTGACAAATTAATGACAACAATTGCGCCAGAAATCCCAGCAGTGCGCCACCATACTAGAAAACGTAGAGCAGGATGAATCCATTGTAGAGTGCACCAATCAGTTTCTTGGACTTTAATGAGTCATAAGCGAGCAAAGATACTGGGACGGCACACTGTTACCGATGACGAGAACACGGCCGCCACCATCAAACCACCCACGTAAGCCGAAACGCGCCGACCAGCCCGGTCAAACGCCGCCAGACGGGCAAAACCAGCTGTCATCATCCGCACTTCGTCAATGCCATTTATCAATCTGCGTGATACGTGATCACGCATGGTCAGTGCCAGTAATTCCCTATGTCCCCTGGGAAACCCAGATTTCGGGATCAGAATCGCCATAATGCCGGAATAGCGGGTTTTCTTCCAGCATCTACAGCGCACCACCAAAACCAGACTCAGTTACCGAGTACGCGATGTTAGCGTGTCTTCACGAGAAGATAAGTTCACTTTCTGACTCCGCTTGCTCGGCCTGCAACTGTGCCATCAACAGCCACGCAGCAATCCTTCACGGATGTTTTGCAGATAAACAGTGATATTTCGAACGATAATACTGCCTTCGCCTGGGTTAAACGTTTTACACCATAACAGACGCTGGGCGGCACCGGTCAGAGAACCCGATAACCATCACCGTTAAAACGTCTTAGCGAACGGTTTGCTTTCAGATTAATCCGCAATGCGCCGTGGAAGGCACCACGCGCAGCTGTTTCGGTTGAAAACAGCTTTGGTTTGCGATTTACCGGTTTGGCGAGGATAACAATTTCTGGCGCGTGTTTTGTAGTGCGTCAGAACCTCAGCCTCCGCTTCTAACAGCCGATTTACGGGAACACCGCGATCAAGACGCGAGCCAAATTGCCGAAAATATAATTGACTTGCATATCGCGCTTATTCATCTAGCACGCCTCCATGCAAGTGACCATCGCTCAGGGTGAGCATGCGATAGGAACGCCGGGAGGTTACGATTGATCATTGTGCGTTGCCATCAATACGGTTACCCAACGCGGTTAAACTCTTCCAAACAGACGTAAAATGCCTTTCTCGACAAGTCCTGAGGTTACCAGTCGGTTCTGTCCGCCAGCAGTAACCCCTTAGTTTCATCCACCGCGCGGGTGGCAATTGCCAACACGCTGTTGTTGGCACCGCCCGAAAGCTGAATAGGGAAGTTGGTTCTCGCTTTGTCCAGCCCGACTTTATCCAGCGCCGCCGACACCCGGCGACCCAATGTCACCGCTGGCAGGCGATAATCAAGCGGCGTAAAAATTGCCACTTAGACAGTAATGTTGGATAATAGATGATCCTGGAAAACCGCTTGCCAATGGGCATAGAAACGAACTTCACGGTTTTCGACCGCGAGCGGTGATGTCATATGAGGACAACCAGATTTTCCCGGCGCTGGGCCGCTCAATCGTACAGGATCAGCTTCGCAGATTTTCTTATTTAGAATACGAATTCCAAAACGCCATCGCCGATGGGTGCATATGGAAATGGTAACGCCCTGCGCAGCGCCTGTACTCCCACCGAGATAATGTGCTGACATGTTCAAAGCGAATCATTGTTAATCCTCTCGGGCAAAATTGCCTCTATAAAGTCGTCCGCCTTAAACGGACGCAAATCCTCTAATACGTTTGCCGACACCAATGTAGCGGATAGGGATACCAAACTGGTCAGCCACCCAGAAAATTACCGACTTCGCTTTGCCGCTTTAAGTTTCGTCAGGCGCGTGGACTCATGTCATAAGCCAACGGCTTCATCGAACAGTGTTGCCTGCTTACCGCGATTCTGCCCGGTGCTGGCATCAATAATTACCAGCATAATTTCATGCGGCGCTTTCAACGTCGAGTTTCTTCATCACGCGGACGATTTTCTTCAACTCTTCTACATCAGGTGCGATTTGTTTCTGCAGGCGCAGTCCGGGGCCTGTAATGGCAATCAGGAAAGTCGATATTGACGCGCTATTTAGCTGCCTGAACCTAGTCGAAAGATAACAGAGGCGAATCTCCGGTATGCTGGGCAATCACCGGAATATTGTTGCGCTGACCCCAGATATGAAGCTGTTCAACGCAGCTTGCACGGAAGATGACCTGCCGCCAGCATCACCGATTTACCCTGCTGCTCAAACTGACGCGCCAGCTTACCAATCGTCGTGGTTTTCCCACACCGTTGAAGCCCACCATCAGCAGTAACCAAACGGCGCTTTGCCTTCAAATATTCAGCGGCTCATCGACTTTCGCCAGAATCTCGCCCATCTCTTCTTCATAGCAAGGCCATAGAGCTGCCTCGGGGTTACGCAACTATGCTGCGCGTACATTTCTGCCTCCGTCAGATTGGTAATTTTACGTGTGGTTTCCACACCCAATTAGGCAATCAAAACAAGCTTAGCTCTTCCAGCTCCTCAAAATACGACCCGGTAGCGTTACGATTTTACCGCGGAACAGGCTGATAAAATCCGAGAACCGAGATTTTCTTTGGTTTTAACAGGCTGCGTTTCAGGCGCGCGAAAATTAATTCTTTGTATTGGTTTTTCCTGCTCCTGAGCGATTTCTTCACCGGCTGCTCTTCTTCTGCCGGAGGAACAACAAGAACACCTCTTCTGCCGCTTCGGCAGCCAGCGGGTTTCCTATGTTCTCGTCGGTGATTTCTTCTTTAGCCGCTTCTTCTTCAGCCGCTATTCGACAATCTCTACGGTTTCGCTTCAGCCTGCCACTCTTCTGAATACGATTTGCTTCGGCGTTGACGTCTTCTAGCGGCAACGGCATTAGCCTCTATTTCTACGTTCAAGCGCTCCGATTTCTTCTACGACAAGCTCCCTGTTGCGTCAAAGGCTACATCTTCCGCTTCAGGCTTCGCTTTTCACTTTCAGCAACCTGTTCAGGTGACTTCTCCACAACGTCGGCAGGGACAAAGTTTCTTCGCTCGGCTTCAGTATGCGCTTGCGGCTGCTCTTCAACCGCTTGTTCAGAGGCCTTCACAGGCTCTATTGCGCCTGAACGATTTCTTCTACAACCGGTTGTTCATATTCTGATTCTGTCTCTTTTCCGGGGTCTGCTCTTTGACCAAGCCCAGCCAAGAAGGCTTTTCTTTTCACATACTGACTTTACAGCCTCCTATGTTGCTTTCATGGCACAGCGTCAAACGCTATGTACATAGCAGCTAAAATGATGAAATAGTCTATCACTTAACTTAATTCACATCAAGCTGCAATATGTTATCTGGCGGATTGAGCAATTTATCATGAAAATGGCAAATCATTCCGGACACAGCGGCCAAATCGCATTATTGGCGGGCATGGCGGGGCTTACTCCTCCCGGTTCCTGTCACTAAGGTCTCTGCCCCACCAATGACTCGGTACGGGGCAACGGTTGTTTAACTGGCTGAGCCGGTTGATGTTGACTCGCCCAATGTCTGGATTGCTTCGCCGGGAGCGGCTCGAGGGTGGCTGAAGGGTTTGGGCGCTACATTGCTCGCGGGGGCAACGGTTGACTTGAGATGGATCGCGTTGGCCATCGAAGTTAATTAAGAATCTGGCGACACTAAAATTAAGAGAGGCAATGCGAAGAGCGTGGTGAACAGCAACGCGATGTCCAATTCCTGGGGGCAAAGGTACACCGCATAATACTGTGTTTGTCGATCACCGTTCGCCTTGGCTTTAGAGACGATAAATTTACTGAGAAGATAGACCGTAGCTGGCTGACGAAGCCCTGTAAGTTGTATGTCGAAAGCGAAGTCGAAAACGGTCTGCCACTGTTCGAGCAAACTGGTCTATTACATCGCAAAATTAGGGTCAGGTAGGCTTATCGGCTGTATCAACGCGACAGCACAAGGAGAAAGTGATGCTGTCGTGATAATATTATGCTTTTGTTAATGCTCTGCGTTTGGGGATTTTAATCCTCAACCTGGGTCGCATCCTTCCCCACGCCCGCTGAATATCTTCTAGTTAACGTGGCGCTGATTTTAAGTTTGGTTTATGG";
    char* pattern = "TTCGATCCAGTACTACGAT";
//    char *text = "Ovo je random dugi string cisto da isprobam kako radi kada imam nesto duze sto mozda malo duze treba parsirati. Parsiranje je jako kul ali ne radi najbolje.";
//    char *pattern = "randikdso";
//    char *text = "banana";
//    char *pattern = "ana";
    int textLength = strlen(text);
    int patternLength = strlen(pattern);
    int k = 8;


    NaiveMismatcher mismatcher(text, textLength, pattern, patternLength, k);
    clock_t begin = clock();
    int* mismatches = mismatcher.findMismatches();
    clock_t end = clock();
    double elapsed_secs1 = double(end - begin) / CLOCKS_PER_SEC;

    for (int i = 0; i < textLength - patternLength + 1; i++) {
        if(0 <= mismatches[i] && mismatches[i] <= k) {
            cout << i << ":" << mismatches[i] << endl;
        }
    }
    vector<StringPointer> subsets = SuffixArray(text, pattern, textLength, patternLength).subsets;
    KangarooMismatcher kangarooMismatcher(subsets, text, textLength, pattern, patternLength, k);
    begin = clock();
    int* kangarooMismatches = kangarooMismatcher.findMismatches();
    end = clock();
    double elapsed_secs2 = double(end - begin) / CLOCKS_PER_SEC;
//
//
    for (int i = 0; i < textLength - patternLength + 1; i++) {
        if(0 <= kangarooMismatches[i] && kangarooMismatches[i] <= k) {
            cout << i << ":" << kangarooMismatches[i] << endl;
        }
    }
    cout << elapsed_secs1;
    cout << "\n";
    cout << elapsed_secs2;
//    cout << "\n";
//    cout << kangarooMismatcher.rmqTime;


    return 0;
}


//int
//main(int argc, char* argv[]) {
//  program_name = argv[0];
//
//  ifstream is;
//
//  // Handle arguments
//  int opt;
//  int option_index;
//  const char* optchar*="hf:s:";
//  while ((opt = getopt_long_only(argc, argv, optchar*, long_options,
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
//    wordchar* word;
//    is >> word;
//    //word.push_back(numeric_limits<int>::max());
//    word.push_back("~~~~~~~~~~~~~");
//    sa_char* tree(word);
//    is.close();
//
//  sa_char*::size_type result;
//    cerr << "STREE: SUFFIXARRAY Ready!" << endl;
//    char input[1000000];
//    while ( cin.getline(input, 1000000) ) {  // read the char* and write it back
//      //cerr << "STree: RECIEVED ***" << input << "***"<< endl;
//    wordchar* totest;
//      ichar*stream a (char*(input), ichar*stream::in);
//      a >> totest;
//      result = tree.find_all_positions_count(totest);
//      //cerr << "SARRAY: Sending result: ***"<<result <<"***"<< endl;
//      cout << result << endl;
//    }
//}
// end of file: main.cpp
