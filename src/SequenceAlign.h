#pragma once
#include <iostream>
#include <vector>
#include "FastaFile.h"
#include <string>
#include <algorithm>
class SequenceAlign {
public:
    SequenceAlign(const FastaFile& fileOne, const FastaFile& fileTwo);

private:

    //main functionality for aligning sequences
    void Align(const std::string& s1, const std::string& s2);

    //output the final result to file
    void Output();

    //output infos
    std::string mHeaderOne;
    std::string mHeaderTwo;
    std::string mFileName;

    //result string
    std::string mResultS1;
    std::string mResultS2;
    int mScore;
};

