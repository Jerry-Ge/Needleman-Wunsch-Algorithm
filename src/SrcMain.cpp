#include "SrcMain.h"
#include <iostream>
#include "FastaFile.h"
#include "DnaTranslator.h"
#include "SequenceAlign.h"
#include <string>
void ProcessCommandArgs(int argc, const char* argv[])
{
	if (argc < 2) {
		std::cerr << "No enough arguments" << std::endl;
	}

	//amino acid counts
	if (argc == 2) {
		FastaFile myFile(argv[1]);

		//initialize the object and perform functions
		DNATranslator translator(myFile);
	}
	//sequence alignment
	else if (argc == 3) {
		FastaFile myFileOne = argv[1];
		FastaFile myFileTwo = argv[2];

		//call the sequence alignment object
		SequenceAlign aligner(myFileOne, myFileTwo);
	}

}
