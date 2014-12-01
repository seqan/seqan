#ifndef APPS_BS_TOOLS_BISAR_BASE_H_
#define APPS_BS_TOOLS_BISAR_BASE_H_

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
    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(fileName)))
		return false;

	String<Dna5Q> seq;
	CharString _id;

    while (!atEnd(seqFileIn))
    {
        readRecord(_id, seq, seqFileIn);
        cropAfterFirst(_id, IsBlank());
        appendRead(store, seq, _id);
    }

    return true;
}


template <typename TFSSpec, typename TFSConfig, typename TFileName>
bool loadReadsCroppedId(FragmentStore<TFSSpec, TFSConfig> & store,
                        TFileName & fileNameL, TFileName & fileNameR)
{
    seqan::SeqFileIn seqFileInL, seqFileInR;
	if (!open(seqFileInL, toCString(fileNameL)) ||
	    !open(seqFileInR, toCString(fileNameR)))
		return false;

	String<Dna5Q> seqL, seqR;
	CharString _idL, _idR;

    while (!atEnd(seqFileInL) && !atEnd(seqFileInR))
    {
        readRecord(_idL, seqL, seqFileInL);
        cropAfterFirst(_idL, IsBlank());
        readRecord(_idR, seqR, seqFileInR);
        cropAfterFirst(_idR, IsBlank());

		appendMatePair(store, seqL, seqR, _idL, _idR);
    }

	SEQAN_ASSERT(atEnd(seqFileInL) && atEnd(seqFileInR));

	return true;
}


#endif
