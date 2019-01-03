#pragma once

#include <string>
#include <fstream>
#include <iostream>

class FastaFile
{
public:
	FastaFile(const char* fileName);

	std::string GetSequence() const;
	std::string GetFileName() const;
	std::string GetHeader() const;

private:
	std::string mHeader;
	std::string mSequence;
	std::string mOutputName;

	void ReadFile(std::ifstream& infile);
};

