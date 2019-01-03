#include "SequenceAlign.h"


/*
 * Constructor
*/
SequenceAlign::SequenceAlign(const FastaFile& fileOne, const FastaFile& fileTwo)
{

    mScore = 0;
    mHeaderOne = fileOne.GetHeader();
    mHeaderTwo = fileTwo.GetHeader();
    mFileName = fileOne.GetFileName();
    //get the actual sequence for comparison
    std::string s1 = fileOne.GetSequence();
    std::string s2 = fileTwo.GetSequence();

    //main calculating function 
    Align(s1, s2);

    //output the result strings to file
    Output();

}


/*
 * main functionality for aligning sequences
*/
void SequenceAlign::Align(const std::string& s1, const std::string& s2)
{
	//create first grid of short of scores
	int m = s2.size() + 1; //row size: size of B + 1
	int n = s1.size() + 1; //col size: size of A + 1

	std::vector<std::vector<short>> scores(m, std::vector<short>(n, 0));

	//create second grid of chars of directions
	std::vector<std::vector<char>> directions(m, std::vector<char>(n, 0));

	//an enum to represent directions
	enum Direction : char
	{
		above,  //0
		aboveLeft, //1
		left //2
	};

	//grid 1 initialization
	scores[0][0] = 0;
	for (int j = 1; j < n; ++j) { //init first row
		scores[0][j] = scores[0][j - 1] - 1;
	}
	for (int i = 1; i < m; ++i) { //init first column
		scores[i][0] = scores[i - 1][0] - 1;
	}

	//grid 2 initialization [arrows]
	for (int j = 1; j < n; ++j) { //init direction for first row
		Direction d = left;
		directions[0][j] = d;
	}
	for (int i = 1; i < m; ++i) { //init direction for first coloumn
		Direction d = above;
		directions[i][0] = d;
	}


	//main dp double loop
	for (int i = 1; i < m; ++i) {
		for (int j = 1; j < n; ++j) {
			char a = s1[j - 1];
			char b = s2[i - 1];

			//three candidate values
			short aboveScore = static_cast<short>(scores[i - 1][j] - 1);
			short aboveLeftScore = (a == b) ? static_cast<short>(scores[i - 1][j - 1] + 1) : static_cast<short>(scores[i - 1][j - 1] - 1);
			short leftScore = static_cast<short>(scores[i][j - 1] - 1);

			//get the max score
			short maxScore = std::max(aboveScore, std::max(aboveLeftScore, leftScore));
			if (maxScore == aboveLeftScore) {
				scores[i][j] = aboveLeftScore;
				Direction d = aboveLeft;
				directions[i][j] = d;
			}
			else if (maxScore == leftScore) {
				scores[i][j] = leftScore;
				Direction d = left;
				directions[i][j] = d;
			}
			else {
				scores[i][j] = aboveScore;
				Direction d = above;
				directions[i][j] = d;
			}
		}
	}

	//get the optimal score
	mScore = scores[m - 1][n - 1];

	//follow the arrows to print out result string
	int i = m - 1, j = n - 1;
	while (i > 0 && j > 0) {
		if (directions[i][j] == aboveLeft) {
			mResultS1 += s1[j-1];
			mResultS2 += s2[i-1];
			j--;
			i--;
		}
		else if (directions[i][j] == left) {
			mResultS1 += s1[j-1];
			mResultS2 += "_";
			j--;
		}
		else { //above case
			mResultS2 += s2[i-1];
			mResultS1 += "_";
			i--;
		}
	}

	//print out remaining chars
	while (j > 0) {
		mResultS1 += s1[j-1];
		mResultS2 += "_";
		j--;
	}
	while (i > 0) {
		mResultS2 += s2[i-1];
		mResultS1 += "_";
		i--;

	}

    //reverse the strings to get the actual solution
    std::reverse(mResultS1.begin(), mResultS1.end());
    std::reverse(mResultS2.begin(), mResultS2.end());
	
}


/*
 * output the final result to file
*/
void SequenceAlign::Output() 
{
	std::string outputName = "match.result";
	std::ofstream ofile(outputName);

	//print out headers and score
    ofile << "A: " << mHeaderOne.substr(1, mHeaderOne.size() - 1) << std::endl;
    ofile << "B: " << mHeaderTwo.substr(1, mHeaderTwo.size() - 1) << std::endl;
    ofile << "Score: " << mScore << std::endl;

    ofile << std::endl;

    int totalLen = mResultS1.size();

	//every line has max of 70 chars
	int level = 0;
	int maxLevel = totalLen / 70;

	//print all full levels
	for (int i = 0; i < maxLevel; ++i) {
		int startIndex = level * 70;
		std::vector<bool> same(70, false);
		for (int i = startIndex; i < startIndex + 70; ++i) {
			if (mResultS1[i] == mResultS2[i]) {
				same[i - startIndex] = true;
			}
			ofile << mResultS1[i];
		}
		ofile << std::endl;

		for (int i = 0; i < 70; ++i) {
			if (same[i]) {
				ofile << "|";
			}
			else {
				ofile << " ";
			}

		}
		ofile << std::endl;
		for (int i = startIndex; i < startIndex + 70; ++i) {
			ofile << mResultS2[i];
		}
		ofile << std::endl;
		ofile << std::endl;

		//increment the index of full level
		level++;
	}

	//print last level
	int startIndex = maxLevel * 70;
	int remainingSize = totalLen - maxLevel * 70;
	std::vector<bool> same(remainingSize, false);
	for (int i = startIndex; i < startIndex + remainingSize; ++i) {
		if (mResultS1[i] == mResultS2[i]) {
			same[i-startIndex] = true;
		}
		ofile << mResultS1[i];
	}
	ofile << std::endl;
	for (int i = 0; i < remainingSize; ++i) {
		if (same[i]) {
			ofile << "|";
		}
		else {
			ofile << " ";
		}
		
	}
	ofile << std::endl;
	for (int i = startIndex; i < startIndex + remainingSize; ++i) {
		ofile << mResultS2[i];
	}
	ofile << std::endl;
}


