#include "FastaFile.h"

FastaFile::FastaFile(const char* fileName)
{
	std::ifstream infile(fileName, std::ios::in | std::ios::ate);

	//pre-allocate space for input file
	size_t fileSize = 0;
	if (infile.is_open())
	{
		fileSize = infile.tellg();
	}
	infile.clear();
	infile.close();
	this->mSequence.reserve(fileSize);

	std::string outputFilePrex(fileName);
	int pos = outputFilePrex.find_first_of('_');
	outputFilePrex = outputFilePrex.substr(6, pos-6);
	this->mOutputName = outputFilePrex;

	//actually read the file
	std::ifstream infile2(fileName);
	ReadFile(infile2);
}


//the actual parser for reading the file
void FastaFile::ReadFile(std::ifstream& infile)
{
	std::getline(infile, this->mHeader); //read the header

	std::string line;
	while (static_cast<bool>(std::getline(infile, line))) {
		this->mSequence.append(line);
		line.clear();
	}
}

std::string FastaFile::GetSequence() const
{
	return this->mSequence;
}

std::string FastaFile::GetFileName() const
{
	return this->mOutputName;
}

std::string FastaFile::GetHeader() const
{
	return this->mHeader;
}