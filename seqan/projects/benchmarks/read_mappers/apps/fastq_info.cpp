/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ==========================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  ==========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ==========================================================================
  A simple analysis tool for sequence files.

  Called with an arbitrary number of FASTA/FASTQ files, it prints
  statistics about the sequence lengths and qualities.

  Usage: fastq_info FILE1.FASTQ [FILE2.FASTQ [FILE3.FASTQ ...]]
  ==========================================================================*/

#include <iostream>
#include <algorithm>
#include <map>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace seqan;


// Computes the median of a sequence of numbers, using the usual
// definition for an even number of elements.
template <typename T>
double median(String<T> & seq) {
    typedef typename Iterator<String<T>, Standard>::Type TIterator;
    TIterator first = begin(seq, Standard());
    TIterator last = end(seq, Standard());

    if (length(seq) % 2 == 1) {
        std::nth_element(first, first + (last - first) / 2, last);
        return *(first + (last - first) / 2);
    } else {
        std::nth_element(first, first + (last - first) / 2, last);
        double result = *(first + (last - first) / 2);
        result += *std::min_element(first + (last - first) / 2, last);
        return result / 2;
    }
}


// Perform as imple analysis of the sequences in the given string set
// with the given ids.  The filename is only given for output
// purposes.
void performAnalysis(char const * filename,
                     StringSet<String<Dna5Q> > const & sequences,
                     StringSet<CharString> const &) {
    typedef StringSet<String<Dna5Q> > TSeqStringSet;
    typedef Iterator<TSeqStringSet const, Standard>::Type TSeqStringSetIterator;
    typedef Iterator<String<Dna5Q> const, Standard >::Type TStringIterator;

    // Compute simple statistics about the sequences.
    String<size_t> sequenceLengths;
    String<int> qualities;
    int minQuality = 1000;
    int maxQuality = 0;
    std::map<int, size_t> qualityHistogram;
    __int64 qualitySum = 0;
    size_t minSequenceLength = ~0;
    size_t maxSequenceLength = 0;
    size_t sequenceLengthSum = 0;
    for (TSeqStringSetIterator it = begin(sequences, Standard()); it != end(sequences, Standard()); ++it) {
        appendValue(sequenceLengths, length(*it));
        sequenceLengthSum += length(*it);
        minSequenceLength = _min(minSequenceLength, length(*it));
        maxSequenceLength = _max(maxSequenceLength, length(*it));

        for (TStringIterator itBasePair = begin(*it, Standard()); itBasePair != end(*it, Standard()); ++itBasePair) {
            qualityHistogram[getQualityValue(*itBasePair)] += 1;
            appendValue(qualities, getQualityValue(*itBasePair));
            qualitySum += getQualityValue(*itBasePair);
            minQuality = _min(minQuality, getQualityValue(*itBasePair));
            maxQuality = _max(maxQuality, getQualityValue(*itBasePair));
        }
    }
    double avgSequenceLength = 1.0 * sequenceLengthSum / length(sequences);
    double avgQuality = 1.0 * qualitySum / sequenceLengthSum;
    double medianQuality = median(qualities);
    double medianSequenceLength = median(sequenceLengths);

    // Print the analysis.
    std::cout << "File: " << filename << std::endl;
    std::cout << "number of sequences: " << length(sequences) << std::endl;
    std::cout << std::endl;

    std::cout << "Sequence Lengths" << std::endl;
    std::cout << "minimal sequence length: " << minSequenceLength << std::endl;
    std::cout << "average sequence length: " << avgSequenceLength << std::endl;
    std::cout << "median  sequence length: " << medianSequenceLength << std::endl;
    std::cout << "maximal sequence length: " << maxSequenceLength << std::endl;
    std::cout << std::endl;

    std::cout << "Base Pair Qualities" << std::endl;
    std::cout << "minimal quality: " << minQuality << std::endl;
    std::cout << "average quality: " << avgQuality << std::endl;
    std::cout << "median  quality: " << medianQuality << std::endl;
    std::cout << "maximal quality: " << maxQuality << std::endl;
    std::cout << std::endl;

    std::cout << "Quality Histogram" << std::endl;
    for (std::map<int, size_t>::const_iterator it = qualityHistogram.begin(); it != qualityHistogram.end(); ++it) {
        printf("%2d\t%8lu\n", it->first, it->second);
    }
    std::cout << std::endl;
    std::cout << std::endl;
}


int main(int argc, char ** argv) {
    // Check arguments.
    if (argc < 2) {
        std::cerr << "Invalid character count!" << std::endl;
        std::cerr << "USAGE: fast_qinfo FILE.FASTQ+" << std::endl;
        return 1;
    }

    // Read in files.
    for (int argi = 1; argi < argc; ++argi) {
        std::cout << "Reading file " << argv[argi] << std::endl;
        MultiSeqFile multiSeqFile;
        if (!open(multiSeqFile.concat, argv[argi], OPEN_RDONLY)) {
            std::cerr << "Could not open file " << argv[argi] << std::endl;
            return 1;
        }
        AutoSeqFormat format;
        guessFormat(multiSeqFile.concat, format);
        split(multiSeqFile, format);
        unsigned seqCount = length(multiSeqFile);
        StringSet<String<Dna5Q> > seqs;
        StringSet<CharString> seqIDs;
        reserve(seqs, seqCount, Exact());
        reserve(seqIDs, seqCount, Exact());
        String<Dna5Q> seq;
        CharString qual;
        CharString id;
        for (unsigned i = 0; i < seqCount; ++i) {
            assignSeq(seq, multiSeqFile[i], format);    // read sequence
            assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
            assignSeqId(id, multiSeqFile[i], format);   // read sequence id
            // Convert ascii to values from 0..62 store dna and quality
            // together in Dna5Q.
            for (unsigned j = 0; j < length(qual) && j < length(seq); ++j)
                assignQualityValue(seq[j], (int)(ordValue(qual[j]) - 33));
            // We use reserve and append, as assign is not supported by
            // StringSet<..., Owner<ConcatDirect<> > >.
            appendValue(seqs, seq, Generous());
            appendValue(seqIDs, id, Generous());
        }

        // Do some analysis.
        performAnalysis(argv[argi], seqs, seqIDs);
    }
    return 0;
}
