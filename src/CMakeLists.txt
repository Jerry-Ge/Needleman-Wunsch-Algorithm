# If you create new headers/cpp files, add them to these list!
set(HEADER_FILES
	SrcMain.h
	FastaFile.h
	DnaTranslator.h
	SequenceAlign.h
)

set(SOURCE_FILES
	SrcMain.cpp
	FastaFile.cpp
	DnaTranslator.cpp
	SequenceAlign.cpp
)

# Don't change this
add_library(src ${SOURCE_FILES} ${HEADER_FILES})
