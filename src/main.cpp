/*
 *
 * String matching with k mismatches
 *
 * Bioinformatics project
 *
 * Authors:
 *  Domagoj Korman
 *  David Romic
 *  Vilim Stubican
 *
 *
 */

#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <math.h>
#include <map>
#include <set>
#include <unordered_map>
#include "sais.h"


using namespace std;


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

        while (j + l < patternLength) {
            int tmpL = 0;
            int tmpI = i;
            int tmpJ = j;
            while (tmpJ < patternLength && text[tmpI++] == pattern[tmpJ++]) {
                tmpL++;
            }

            if (tmpL > l) {
                l = tmpL;
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
    int* suffixArray;
    int* inverseSuffixArray;
    int* lcp;

    int updateMismatches(int positionInText);

    void init();

    void createLCP();

    int RMQ(int index1, int index2);

public:
    KangarooMismatcher(int sa[], int isa[], char *text, int textLength, char *pattern, int patternLength, int kVal);
};

KangarooMismatcher::KangarooMismatcher(int sa[], int isa[], char *text, int textLength, char *pattern, int patternLength,
                                       int kVal)
        : Mismatcher(text, textLength, pattern, patternLength, kVal) {
    this->suffixArray = sa;
    this->inverseSuffixArray = isa;
    init();
}

void KangarooMismatcher::init() {
    cout << "\nInit Kangaroo mismatcher\n";
    createLCP();
    cout << "\nCreated LCP, algorithm starting\n";
}

void KangarooMismatcher::createLCP() {
    char* textAndPattern = new char[textLength + patternLength + 1];
    strcpy(textAndPattern, text);
    strcat(textAndPattern, "~");
    strcat(textAndPattern, pattern);
    lcp = new int[textLength + patternLength + 1];
    for (int i = 0; i < textLength + patternLength; ++i) {

        //Calculate LCP between string and string i + 1
        int pos1 = suffixArray[i];
        int pos2 = suffixArray[i + 1];
        int length1 = patternLength + textLength + 1 - pos1;
        int length2 = patternLength + textLength + 1 - pos2;
        int lcpValue = 0;
        int j = 0;
        while (j < length1 && j < length2
               && textAndPattern[pos1 + j] == textAndPattern[pos2 + j]) {
            lcpValue++;
            j++;
        }
        lcp[i] = lcpValue;
    }
    lcp[textLength + patternLength] = 0;
}

int KangarooMismatcher::RMQ(int index1, int index2) {
    int lcpIndex1 = inverseSuffixArray[index1];
    int lcpIndex2 = inverseSuffixArray[index2];
    int rmq = 100000;
    if (lcpIndex1 > lcpIndex2) {
        for (int i = lcpIndex2; i < lcpIndex1; ++i) {
            if (lcp[i] < rmq) {
                rmq = lcp[i];
                if (rmq == 0) return 0;
            }
        }
        return rmq;
    }
    if (lcpIndex2 > lcpIndex1) {
        for (int i = lcpIndex1; i < lcpIndex2; ++i) {
            if (lcp[i] < rmq) {
                rmq = lcp[i];
                if (rmq == 0) return 0;
            }
        }
        return rmq;
    }
    return -2000;
}

int KangarooMismatcher::updateMismatches(int positionInText) {
    if (positionInText <= textLength - patternLength) {
        int mismatches = 0;
        int l = 0;

        while (true) {
            l += RMQ(positionInText + l, textLength + l + 1) + 1;
            if (l > patternLength) {
                return mismatches;
            }
            if (l == patternLength) {
                return mismatches + 1;
            }
            mismatches++;
            if (mismatches > k) {
                return mismatches;
            }
        }
    } else {
        return k + 1;
    }
}

int main(void) {

    //CHANGE TEXT AND PATTERN HERE
    char text[] = "banana";
    char pattern[] = "ana";
    int k = 1;


    int text_len = strlen(text);
    int patternLength = strlen(pattern);
    char* textAndPattern = new char[text_len + patternLength + 1];
    strcpy(textAndPattern, text);
    strcat(textAndPattern, "~");
    strcat(textAndPattern, pattern);
    int *suffix_array;
    int *inverse_sa;

    suffix_array = (int *) malloc((text_len + patternLength + 1) * sizeof(int));
    inverse_sa = (int *) malloc((text_len + patternLength + 1) * sizeof(int));

    if (sais((unsigned char *) textAndPattern, suffix_array, text_len + patternLength + 1) != 0) {
        printf("sais failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i=0; i < text_len + patternLength + 1; i++) {
        inverse_sa[suffix_array[i]] = i;
    }

    NaiveMismatcher mismatcher(text, text_len, pattern, patternLength, k);
    clock_t begin = clock();
    int* mismatches = mismatcher.findMismatches();
    clock_t end = clock();
    double elapsed_secs1 = double(end - begin) / CLOCKS_PER_SEC;

    for (int i = 0; i < text_len - patternLength + 1; i++) {
        if(0 <= mismatches[i] && mismatches[i] <= k) {
            cout << i << ":" << mismatches[i] << endl;
        }
    }

    begin = clock();
    KangarooMismatcher kangarooMismatcher(suffix_array, inverse_sa, text, text_len, pattern, patternLength, k);
    int* kangarooMismatches = kangarooMismatcher.findMismatches();
    end = clock();
    double elapsed_secs2 = double(end - begin) / CLOCKS_PER_SEC;

    for (int i = 0; i < text_len - patternLength + 1; i++) {
        if(0 <= kangarooMismatches[i] && kangarooMismatches[i] <= k) {
            cout << i << ":" << kangarooMismatches[i] << endl;
        }
    }

    cout << endl << "Time:" << endl;
    cout << elapsed_secs1 << endl;
    cout << elapsed_secs2 << endl;
    return 0;
}


