#ifndef CORE_APPS_BS_TOOLS_BISAR_BASE_H_
#define CORE_APPS_BS_TOOLS_BISAR_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fstream>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/bam_io.h>
#include <seqan/score.h> 

using namespace std;
using namespace seqan;


template <typename TFSSpec, typename TFSConfig, typename TFileName>
bool loadReadsCroppedId(FragmentStore<TFSSpec, TFSConfig> &store, TFileName &fileName)
{
	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY))
		return false;

	// guess file format and split into sequence fractions
	AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

	// reserve space in fragment store
	unsigned seqOfs = length(store.readStore);
	unsigned seqCount = length(multiSeqFile);
	reserve(store.readStore, seqOfs + seqCount);
	reserve(store.readSeqStore, seqOfs + seqCount);
	reserve(store.readNameStore, seqOfs + seqCount);

	// read sequences
	String<Dna5Q> seq;
	CharString qual;
	CharString _id;

	for (unsigned i = 0; i < seqCount; ++i)
	{
		assignSeq(seq, multiSeqFile[i], format);    // read sequence
		assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignCroppedSeqId(_id, multiSeqFile[i], format);  // read sequence id up to the first whitespace
		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		// TODO: support different ASCII represenations of quality values
		assignQualities(seq, qual);
		appendRead(store, seq, _id);
	}
    return true;
}


template <typename TFSSpec, typename TFSConfig, typename TFileName>
bool loadReadsCroppedId(FragmentStore<TFSSpec, TFSConfig> & store, TFileName & fileNameL, TFileName & fileNameR)
{
	MultiSeqFile multiSeqFileL, multiSeqFileR;
	if (!open(multiSeqFileL.concat, toCString(fileNameL), OPEN_RDONLY))
		return false;
	if (!open(multiSeqFileR.concat, toCString(fileNameR), OPEN_RDONLY))
		return false;

	// Guess file format and split into sequence fractions
	AutoSeqFormat formatL, formatR;
	guessFormat(multiSeqFileL.concat, formatL);
	split(multiSeqFileL, formatL);
	guessFormat(multiSeqFileR.concat, formatR);
	split(multiSeqFileR, formatR);

    // Check that both files have the same number of reads
	SEQAN_ASSERT_EQ(length(multiSeqFileL), length(multiSeqFileR));

	// Reserve space in fragment store
	unsigned seqOfs = length(store.readStore);
	unsigned seqCountL = length(multiSeqFileL);
	unsigned seqCountR = length(multiSeqFileR);
	reserve(store.readStore, seqOfs + seqCountL + seqCountR);
	reserve(store.readSeqStore, seqOfs + seqCountL + seqCountR);
	reserve(store.readNameStore, seqOfs + seqCountL + seqCountR);

	// Read in sequences
	String<Dna5Q> seq[2];
	CharString qual[2];
	CharString _id[2];

	for (unsigned i = 0; i < seqCountL; ++i) {
		assignSeq(seq[0], multiSeqFileL[i], formatL);    // read sequence
		assignQual(qual[0], multiSeqFileL[i], formatL);  // read ascii quality values
        assignCroppedSeqId(_id[0], multiSeqFileL[i], formatL);  // read sequence id up to the first whitespace 
		assignSeq(seq[1], multiSeqFileR[i], formatR);    // read sequence
		assignQual(qual[1], multiSeqFileR[i], formatR);  // read ascii quality values
		assignCroppedSeqId(_id[1], multiSeqFileR[i], formatR);  // read sequence id up to the first whitespace

		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		// TODO: support different ASCII represenations of quality values
		for (int j = 0; j < 2; ++j)
			assignQualities(seq[j], qual[j]);
		
		appendMatePair(store, seq[0], seq[1], _id[0], _id[1]);
	}
	return true;
}


#endif
