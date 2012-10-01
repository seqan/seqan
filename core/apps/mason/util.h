// ==========================================================================
//                          Mason - A Read Simulator
// ==========================================================================
// Copyright (C) 2010 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Utility Code.
// ==========================================================================

#ifndef UTIL_H_
#define UTIL_H_

using namespace seqan;

// Helper function to write out a {Dna,Dna5}Q sequence with qualities to a
// FASTQ file.
template <typename TStream, typename TIdStringSet, typename TSeqStringSet>
void write(TStream & stream,
           TIdStringSet & seqIds,
           TSeqStringSet & sequences,
           Fastq const &) {
//IOREV definitely move to fastq format module, once we have it
    typedef TSeqStringSet TStringSet;
    typedef typename Position<TStringSet>::Type TPosition;

    CharString qualBuffer;
    for (TPosition i = 0; i < length(sequences); ++i) {
        stream << "@" << seqIds[i] << std::endl;
        stream << sequences[i] << std::endl;
        stream << "+" << /*seqIds[i] << */std::endl;
        resize(qualBuffer, length(sequences[i]), Exact());
        for (TPosition j = 0; j < length(sequences[i]); ++j)
            qualBuffer[j] = getQualityValue(sequences[i][j]) + '!';
        stream << qualBuffer << std::endl;
    }
}

template <typename TStream, typename TIdString, typename TIdStringSpec, typename TSeqString, typename TSeqStringSpec>
void write(TStream & stream,
           StringSet<TIdString, TIdStringSpec> & seqIds,
           StringSet<TSeqString, TSeqStringSpec> & sequences,
           Fasta const &) {
//IOREV definitely move to fastq format module, once we have it
    typedef StringSet<TSeqString, TSeqStringSpec> TStringSet;
    typedef typename Position<TStringSet>::Type TPosition;

    CharString qualBuffer;
    for (TPosition i = 0; i < length(sequences); ++i) {
        stream << ">" << seqIds[i] << std::endl;
        stream << sequences[i] << std::endl;
    }
}

#endif  // UTIL_H_
