#pragma once
#include "FastaFile.h"
#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>

class DNATranslator
{
public:
	DNATranslator(const FastaFile& file);

private:
	void initMap();

	void Translate(const std::string& inputSequence);

	void Print(const std::string& ouputFileName);

	unsigned mTotalCount;
	std::string mHeader;
	std::vector<int> mCounts;
	std::unordered_map<char, std::string> mTranslateName;
	std::unordered_map<std::string, std::pair<char, int> > mToCounts;
	std::unordered_map<int, std::string> mToName;
};
