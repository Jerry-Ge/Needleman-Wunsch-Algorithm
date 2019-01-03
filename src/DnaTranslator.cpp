#include "DnaTranslator.h"

DNATranslator::DNATranslator(const FastaFile& file)
{
	mTotalCount = 0;
	mHeader = file.GetHeader();
	mHeader = mHeader.substr(1, mHeader.size() - 1);
	initMap();

	//allocate space for map
	for (int i = 0; i < 20; ++i) {
		mCounts.push_back(0);
	}
	std::string inputSequence = file.GetSequence();
	std::string outputFilePrefix = file.GetFileName();

	Translate(inputSequence);
	Print(outputFilePrefix);
}

/*
 * The main translation function 
*/
void DNATranslator::Translate(const std::string& inputSequence)
{
	std::vector<std::vector<int>> states{
		//St| T | C | A | G |
		//-------------------
		//0 | 0 | 0 | 1 | 0 | //looking for start codon "A"
			{ 0,  0,  1,  0 },
        //1 | 2 | 0 | 0 | 0 | //looking for start codon "T"
			{ 2,  0,  1,  0 },
		//2 | 0 | 0 | 0 | M | //M and start
			{ 0,  0,  1,  'M'},
		//3 | 4 | 0 | 1 | M | // T**, C**, A**, G**
			{ 4,  9,  14, 19 },
		//4 | TT* | TC* | TA* | TG* | // T**
			{ 5,  6,   7,  8 }, 
		//5 | TTT | TTC | TTA | TTG | // TT*
			{ 'F',  'F', 'L',  'L' },
		//6 | TCT | TCC | TCA | TCG | // TC*
			{ 'S',  'S', 'S',  'S' },
		//7 | TAT | TAC | TAA | TAG | // TA*
			{ 'Y', 'Y',    0,    0 }, 
		//8 | TGT | TGC | TGA | TGG | // TG*
			{ 'C', 'C' ,   0,   'W'}, 
		//9 | CT* | CC* | CA* | CG* | //C**
			{ 10,   11,   12,   13 }, 
		//10| CTT | CTC | CTA | CTG | //CT* 
			{ 'L', 'L',   'L',  'L'}, 
		//11| CCT | CCC | CCA | CCG | //CC* 
			{ 'P',  'P',  'P',  'P', },
		//12| CAT | CAC | CAA | CAG | //CA*
			{ 'H',  'H',  'Q',  'Q' },
		//13| CGT | CGC | CGA | CGG | //CG* 
			{ 'R',  'R',  'R',  'R' }, 
		//14| AT* | AC* | AA* | AG* | //A**
			{ 15 ,  16,   17,   18  },
		//15| ATT | ATC | ATA | ATG | //AT*
			{ 'I',  'I',  'I', 'M'  }, 
		//16| ACT | ACC | ACA | ACG | //AC*
			{ 'T',  'T',  'T',  'T' }, 
		//17| AAT | AAC | AAA | AAG | //AA* 
			{ 'N',  'N',  'K',  'K' },
		//18| AGT | AGC | AGA | AGG | //AG* 
			{ 'S',  'S',  'R',  'R' },
		//19| GT* | GC* | GA* | GG* | //G**
			{ 20,   21,   22,   23  },
		//20| GTT | GTC | GTA | GTG | //GT*
			{ 'V',  'V',  'V',  'V' }, 
		//21| GCT | GCC | GCA | GCG | //GC*
			{ 'A',  'A',  'A',  'A'},
		//22| GAT | GAC | GAA | GAG | //GA* 
			{ 'D',  'D',  'E',  'E'},
		//23| GGT | GGC | GGA | GGG | //GG* 
			{ 'G',  'G',  'G',  'G'}
	};
	


	int inputSize = inputSequence.size();
	int currState = 0;
	for (unsigned i = 0; i < inputSize; ++i) {
		std::vector<int> currVec = states[currState];
		
		//get the index based on current char
		int currIndex;
		if (inputSequence[i] == 'T') {
			currIndex = 0;
		}
		else if (inputSequence[i] == 'C') {
			currIndex = 1;
		}
		else if (inputSequence[i] == 'A') {
			currIndex = 2;
		}
		else {
			currIndex = 3;
		}


		if (currVec[currIndex] < 65) //not a alphabet
		{ 
			currState = currVec[currIndex];
		}
		else { //it's a alphabete
			std::string longName = mTranslateName[currVec[currIndex]];

			int countIndex = mToCounts[longName].second;
			mCounts[countIndex]++;
			mTotalCount++;
			currState = 3;
		}
	}
}

//helper printer function to file
void DNATranslator::Print(const std::string& ouputFileName)
{
	std::string ofileName = "amino.txt";
	std::ofstream ofile(ofileName);
	
	ofile << mHeader << std::endl;

	ofile << "Total amino acids produced: " << mTotalCount << std::endl;

	for (unsigned i = 0; i < mCounts.size(); ++i) {
		std::string longName = mToName[i];
		char shortName = mToCounts[longName].first;
		ofile << "(" << shortName << ") " << longName << ": " << mCounts[i] << std::endl;
	}
}

//helper function to init all map pairs
void DNATranslator::initMap()
{
	mTranslateName['A'] = "Alanine";
	mTranslateName['C'] = "Cysteine";
	mTranslateName['D'] = "Aspartic acid";
	mTranslateName['E'] = "Glutamic acid";
	mTranslateName['F'] = "Phenylalanine";
	mTranslateName['G'] = "Glycine";
	mTranslateName['H'] = "Histidine";
	mTranslateName['I'] = "Isoleucine";
	mTranslateName['K'] = "Lysine";
	mTranslateName['L'] = "Leucine";
	mTranslateName['M'] = "Methionine";
	mTranslateName['N'] = "Asparagine";
	mTranslateName['P'] = "Proline";
	mTranslateName['Q'] = "Glutamine";
	mTranslateName['R'] = "Arginine";
	mTranslateName['S'] = "Serine";
	mTranslateName['T'] = "Threonine";
	mTranslateName['V'] = "Valine";
	mTranslateName['W'] = "Tryptophan";
	mTranslateName['Y'] = "Tyrosine";

	mToCounts["Alanine"] = { 'A', 0 };
	mToCounts["Cysteine"] = { 'C', 1 };
	mToCounts["Aspartic acid"] = { 'D', 2 };
	mToCounts["Glutamic acid"] = { 'E', 3 };
	mToCounts["Phenylalanine"] = { 'F', 4 };
	mToCounts["Glycine"] = { 'G', 5 };
	mToCounts["Histidine"] = { 'H', 6 };
	mToCounts["Isoleucine"] = { 'I', 7 };
	mToCounts["Lysine"] = { 'K', 8 };
	mToCounts["Leucine"] = { 'L', 9 };
	mToCounts["Methionine"] = { 'M', 10 };
	mToCounts["Asparagine"] = { 'N', 11 };
	mToCounts["Proline"] = { 'P', 12 };
	mToCounts["Glutamine"] = { 'Q', 13 };
	mToCounts["Arginine"] = { 'R', 14 };
	mToCounts["Serine"] = { 'S', 15 };
	mToCounts["Threonine"] = { 'T', 16 };
	mToCounts["Valine"] = { 'V', 17 };
	mToCounts["Tryptophan"] = { 'W', 18 };
	mToCounts["Tyrosine"] = { 'Y', 19 };

	mToName[0] = "Alanine";
	mToName[1] = "Cysteine";
	mToName[2] = "Aspartic acid";
	mToName[3] = "Glutamic acid";
	mToName[4] = "Phenylalanine";
	mToName[5] = "Glycine";
	mToName[6] = "Histidine";
	mToName[7] = "Isoleucine";
	mToName[8] = "Lysine";
	mToName[9] = "Leucine";
	mToName[10] = "Methionine";
	mToName[11] = "Asparagine";
	mToName[12] = "Proline";
	mToName[13] = "Glutamine";
	mToName[14] = "Arginine";
	mToName[15] = "Serine";
	mToName[16] = "Threonine";
	mToName[17] = "Valine";
	mToName[18] = "Tryptophan";
	mToName[19] = "Tyrosine";
}