// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Authors: Gianvito Urgese <gianvito.urgese@polito.it>
//          Lily Shellhammer <lily.shellhammer@gmail.com>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_RNA_IO_EBPSEQ_READ_WRITE_H_
#define SEQAN_INCLUDE_SEQAN_RNA_IO_EBPSEQ_READ_WRITE_H_

#include <seqan/sequence.h>
#include <regex>
#include <string>
#include <cstdint>

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Ebpseq
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Ebpseq
 * @headerfile <seqan/rna_io.h>
 * @brief Extended bpseq format for RNA structures (*.ebpseq).
 * @signature typedef Tag<Ebpseq_> Ebpseq;
 * @see FileFormats#RnaStruct
 */
struct Ebpseq_;
typedef Tag<Ebpseq_> Ebpseq;

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Ebpseq, T> :
    public MagicHeader<Nothing, T> {};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Ebpseq, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Ebpseq, T>::VALUE[1] =
{
    ".ebpseq"     // default output extension
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                         [EbpseqHeader]
// ----------------------------------------------------------------------------

template <typename TVec, typename TElem>
inline typename Size<TVec>::Type _findIndex(TVec const & set, TElem const & elem)
{
    typename Size<TVec>::Type idx = 0;
    while (idx < length(set) && set[idx] != elem)
    {
        ++idx;
    }
    return idx;
}

template <typename TForwardIter>
inline void
readHeader(RnaHeader & header, RnaIOContext & context, TForwardIter & iter, Ebpseq const & /*tag*/)
{
    CharString buffer;
    char key;
    clear(header);
    clear(context);

    skipUntil(iter, NotFunctor<IsWhitespace>());
    while (value(iter) == '#' && value(++iter) == '#')
    {
        // read semantics key
        skipOne(iter);
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readOne(key, iter, IsAlpha());

        // read value of header line
        CharString identifier;
        readUntil(identifier, iter, OrFunctor<IsBlank, EqualsChar<':'> >());
        skipOne(iter);
        skipUntil(iter, NotFunctor<IsBlank>());
        readUntil(buffer, iter, IsNewline());

        switch (key)
        {
            case 'G': append(header.description, buffer);
                      break;

            case 'S': appendValue(header.seqLabels, buffer);
                      appendValue(context.seqLabels, buffer);
                      appendValue(context.seqIdent, identifier);
                      break;

            case 'F': appendValue(header.fixLabels, buffer);
                      appendValue(context.fixLabels, buffer);
                      appendValue(context.fixIdent, identifier);
                      break;

            case 'M': appendValue(header.bppLabels, buffer);
                      appendValue(context.bppLabels, buffer);
                      appendValue(context.bppIdent, identifier);
                      break;

            case 'T': appendValue(header.typeLabels, buffer);
                      appendValue(context.typIdent, identifier);
                      break;

            default: SEQAN_THROW(ParseError("ERROR: Unknown key in header."));
        }

        clear(buffer);
        skipUntil(iter, NotFunctor<IsWhitespace>()); // move to next ##
    }
}

// ----------------------------------------------------------------------------
// Function readRecord()                                         [EbpseqRecord]
// ----------------------------------------------------------------------------
// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, RnaIOContext & context, TForwardIter & iter, Ebpseq const & /*tag*/)
{
    if (length(context.seqIdent) >= UINT32_MAX)
    {
        SEQAN_THROW(ParseError("ERROR: Capacity exceeded. Cannot read more than " +
                               std::to_string(UINT32_MAX - 1u) + " sequences."));
    }

    std::string buffer;
    typename Size<StringSet<CharString> >::Type idx;
    clear(record);

    // read sequence and type indices
    skipUntil(iter, NotFunctor<EqualsChar<'#'> >());
    while (value(iter) != '\n' && value(iter) != '\r')
    {
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(buffer, iter, IsWhitespace());
        if (startsWith(buffer, "S"))
        {
            idx = _findIndex(context.seqIdent, buffer.substr(1));
            if (idx < length(context.seqIdent))
            {
                record.name = context.seqLabels[idx];
                record.recordID = static_cast<std::uint32_t>(idx);
            }
            else
            {
                SEQAN_THROW(ParseError("ERROR: Unknown tag symbol " + buffer));
            }
        }
        else if (startsWith(buffer, "T"))
        {
            idx = _findIndex(context.typIdent, buffer.substr(1));
            if (idx < length(context.typIdent))
            {
                appendValue(record.typeID, idx);
            }
            else
            {
                SEQAN_THROW(ParseError("ERROR: Unknown tag symbol " + buffer));
            }
        }
        else
        {
            SEQAN_THROW(ParseError("ERROR: Unknown tag symbol " + buffer));
        }
        clear(buffer);
        skipUntil(iter, NotFunctor<IsBlank>());
    }
    if (record.hasUndefinedID())
    {
        SEQAN_THROW(ParseError("ERROR: Record must contain a sequence tag."));
    }

    // read column semantics
    StringSet<std::string> columnLabels;
    String<unsigned> fixIDs;
    String<unsigned> bppIDs;
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readOne(buffer, iter, EqualsChar<'#'>());
    clear(buffer);
    while (value(iter) != '\n' && value(iter) != '\r')
    {
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(buffer, iter, IsWhitespace());
        appendValue(columnLabels, buffer);
        if (startsWith(buffer, "F"))
        {
            idx = _findIndex(context.fixIdent, buffer.substr(1));
            appendValue(fixIDs, idx);
            if (idx == length(context.fixIdent))
                SEQAN_THROW(ParseError("ERROR: Unknown label " + buffer));
        }
        else if (startsWith(buffer, "M"))
        {
            idx = _findIndex(context.bppIdent, buffer.substr(1));
            appendValue(bppIDs, idx);
            if (idx == length(context.bppIdent))
                SEQAN_THROW(ParseError("ERROR: Unknown label " + buffer));
        }
        clear(buffer);
        skipUntil(iter, NotFunctor<IsBlank>());
    }
    skipUntil(iter, NotFunctor<IsWhitespace>());

    // allocations
    unsigned line {1};
    resize(record.reactivity, length(record.typeID));
    resize(record.reactError, length(record.typeID));
    resize(record.fixedGraphs, length(fixIDs));
    resize(record.bppMatrGraphs, length(bppIDs));

    // set graph specs
    typename Size<String<unsigned> >::Type ii;
    for (ii = 0; ii < length(fixIDs); ++ii)
        record.fixedGraphs[ii].specs = context.fixLabels[fixIDs[ii]];
    for (ii = 0; ii < length(bppIDs); ++ii)
        record.bppMatrGraphs[ii].specs = context.bppLabels[bppIDs[ii]];

    // read matrix
    while (!atEnd(iter) && value(iter) != '#')
    {
        unsigned fixID = 0;
        unsigned bppID = 0;
        forEach(columnLabels, [&](std::string const & label)
        {
            unsigned pairPos;
            typename Size<String<std::size_t> >::Type tID;
            readUntil(buffer, iter, IsWhitespace());

            if (label == "I")
            {
                if (line == 1u && !lexicalCast(record.offset, buffer))
                {
                    SEQAN_THROW(BadLexicalCast(record.offset, buffer));
                }
            }
            else if (label == "NT")
            {
                appendValue(record.sequence, buffer[0]);
            }
            else if (label == "QU")
            {
                appendValue(record.quality, buffer[0]);
            }
            else if (startsWith(label, "F"))
            {
                RnaStructureGraph & graph = record.fixedGraphs[fixID++];
                addVertex(graph.inter);
                if (!lexicalCast(pairPos, buffer))
                {
                    SEQAN_THROW(BadLexicalCast(pairPos, buffer));
                }
                if (pairPos != 0 && line + record.offset - 1 > pairPos)    // add edge if base is connected
                {
                    addEdge(graph.inter, pairPos - record.offset, line - 1, 1.0);
                }
            }
            else if (startsWith(label, "M"))
            {
                RnaStructureGraph & graph = record.bppMatrGraphs[bppID++];
                addVertex(graph.inter);
                if (buffer.find('<') == std::string::npos)
                {
                    SEQAN_THROW(ParseError("ERROR: Expected a bracket pair: < >"));
                }
                if (buffer.find('>') == std::string::npos)
                {
                    readUntil(buffer, iter, EqualsChar<'>'>());
                    skipOne(iter);
                }
                std::regex conn_regex("\\d+/\\d\\.?\\d*");
                std::smatch match;
                while (std::regex_search(buffer, match, conn_regex))
                {
                    double prob;
                    std::string::size_type pos = match.str().find('/');
                    if (!lexicalCast(pairPos, match.str().substr(0, pos)))
                    {
                        SEQAN_THROW(BadLexicalCast(pairPos, match.str().substr(0, pos)));
                    }
                    if (!lexicalCast(prob, match.str().substr(pos + 1)))
                    {
                        SEQAN_THROW(BadLexicalCast(prob, match.str().substr(pos + 1)));
                    }
                    if (pairPos != 0 && line + record.offset - 1 > pairPos)
                    {
                        addEdge(graph.inter, pairPos - record.offset, line - 1, prob);
                    }
                    buffer = match.suffix();
                }
            }
            else if (startsWith(label, "RE"))
            {
                idx = _findIndex(context.typIdent, label.substr(2));
                tID = _findIndex(record.typeID, idx);
                if (tID < length(record.typeID))
                {
                    float react;
                    if (!lexicalCast(react, buffer))
                    {
                        SEQAN_THROW(BadLexicalCast(react, buffer));
                    }
                    appendValue(record.reactError[tID], react);
                }
                else
                {
                    SEQAN_THROW(ParseError("ERROR: Unknown label " + label));
                }
            }
            else if (startsWith(label, "R"))
            {
                idx = _findIndex(context.typIdent, label.substr(1));
                tID = _findIndex(record.typeID, idx);
                if (tID < length(record.typeID))
                {
                    float react;
                    if (!lexicalCast(react, buffer))
                    {
                        SEQAN_THROW(BadLexicalCast(react, buffer));
                    }
                    appendValue(record.reactivity[tID], react);
                }
                else
                {
                    SEQAN_THROW(ParseError("ERROR: Unknown label " + label));
                }
            }
            else
            {
                SEQAN_THROW(ParseError("ERROR: Unknown label " + label));
            }
            skipUntil(iter, NotFunctor<IsBlank>());
            clear(buffer);
        });
        ++line;
        skipLine(iter);
    }
    record.seqLen = length(record.sequence);
}

// ----------------------------------------------------------------------------
// Function writeHeader()                                        [EbpseqHeader]
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeHeader(TTarget & target, RnaHeader const & header, RnaIOContext & context, Ebpseq const & /*tag*/)
{
    typename Size<StringSet<CharString> >::Type idx;
    clear(context);
    if (!empty(header.description))
    {
        write(target, "## G: ");
        write(target, header.description);
        writeValue(target, '\n');
    }

    for (idx = 1; idx <= length(header.seqLabels); ++idx)
    {
        write(target, "## S");
        appendNumber(target, idx);
        write(target, ": ");
        write(target, header.seqLabels[idx - 1]);
        writeValue(target, '\n');
        appendValue(context.seqIdent, std::to_string(idx));
    }

    for (idx = 1; idx <= length(header.fixLabels); ++idx)
    {
        write(target, "## F");
        appendNumber(target, idx);
        write(target, ": ");
        write(target, header.fixLabels[idx - 1]);
        writeValue(target, '\n');
        appendValue(context.fixIdent, std::to_string(idx));
        appendValue(context.fixLabels, header.fixLabels[idx - 1]);
    }

    for (idx = 1; idx <= length(header.bppLabels); ++idx)
    {
        write(target, "## M");
        appendNumber(target, idx);
        write(target, ": ");
        write(target, header.bppLabels[idx - 1]);
        writeValue(target, '\n');
        appendValue(context.bppIdent, std::to_string(idx));
        appendValue(context.bppLabels, header.bppLabels[idx - 1]);
    }

    for (idx = 1; idx <= length(header.typeLabels); ++idx)
    {
        write(target, "## T");
        appendNumber(target, idx);
        write(target, ": ");
        write(target, header.typeLabels[idx - 1]);
        writeValue(target, '\n');
        appendValue(context.typIdent, std::to_string(idx));
    }
}

// ----------------------------------------------------------------------------
// Function createPseudoHeader()
// ----------------------------------------------------------------------------

inline void createPseudoHeader(RnaHeader & header, std::vector<RnaRecord> & records)
{
    if (length(records) >= UINT32_MAX)
    {
        SEQAN_THROW(ParseError("ERROR: Capacity exceeded. Cannot read more than " +
                               std::to_string(UINT32_MAX - 1u) + " sequences."));
    }
    typedef typename Size<std::vector<RnaRecord> >::Type TSizeRnaRecordVector;
    bool hasFixed = false;
    bool hasBppMatr = false;

    for (TSizeRnaRecordVector idx = 0; idx < length(records); ++idx)
    {
        records[idx].recordID = static_cast<std::uint32_t>(idx);
        if (!empty(records[idx].comment))
        {
            if (!empty(header.description))
            {
                append(header.description, "\t");
            }
            append(header.description, records[idx].comment);
        }

        appendValue(header.seqLabels, records[idx].name);

        for (TSizeRnaRecordVector gr = 0; gr < length(records[idx].fixedGraphs); ++gr)
        {
            records[idx].fixedGraphs[gr].specs = "n/a";
            hasFixed = true;
        }

        for (TSizeRnaRecordVector gr = 0; gr < length(records[idx].bppMatrGraphs); ++gr)
        {
            records[idx].bppMatrGraphs[gr].specs = "n/a";
            hasBppMatr = true;
        }
    }

    if (hasFixed)
        appendValue(header.fixLabels, "n/a");

    if (hasBppMatr)
        appendValue(header.bppLabels, "n/a");
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                        [EbpseqRecord]
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, RnaIOContext & context, Ebpseq const & /*tag*/)
{
    if (empty(record.sequence) && length(rows(record.align)) != 1)
    {
        SEQAN_THROW(ParseError("ERROR: Ebpseq formatted file cannot contain an alignment."));
    }
    if (empty(context.seqIdent))       // origin is ebpseq, but header was not printed
    {
        SEQAN_THROW(ParseError("ERROR: Print ebpseq header first."));
    }
    else if (record.recordID >= length(context.seqIdent))
    {
        SEQAN_THROW(ParseError("ERROR: Record ID exceeds number of sequences in ebpseq header."));
    }

    // write head line with S and T tags
    write(target, "# S");
    write(target, context.seqIdent[record.recordID]);
    forEach(record.typeID, [&](std::size_t tID)
    {
        if (tID >= length(context.typIdent))
        {
            SEQAN_THROW(ParseError("ERROR: Type ID exceeds number of types in ebpseq header."));
        }
        write(target, " T");
        write(target, context.typIdent[tID]);
    });
    writeValue(target, '\n');

    // error checks
    if (length(record.fixedGraphs) > length(context.fixIdent))
    {
        SEQAN_THROW(ParseError("ERROR: Graph ID exceeds number of fixed graphs in ebpseq header."));
    }
    if (length(record.bppMatrGraphs) > length(context.bppIdent))
    {
        SEQAN_THROW(ParseError("ERROR: Graph ID exceeds number of bpp matrix graphs in ebpseq header."));
    }

    // write column header
    write(target, "# I\tNT");
    if (!empty(record.quality))
    {
        write(target, "\tQU");
    }
    forEach(record.typeID, [&](std::size_t tID)
    {
        write(target, "\tR");
        write(target, context.typIdent[tID]);

        write(target, "\tRE");
        write(target, context.typIdent[tID]);
    });

    typename Size<StringSet<CharString> >::Type graphIdx;
    forEach(record.fixedGraphs, [&](RnaStructureGraph const & graph)
    {
        graphIdx = _findIndex(context.fixLabels, graph.specs);
        if (graphIdx == length(context.fixLabels))
            SEQAN_THROW(ParseError("ERROR: Fixed structure description not found in header."));

        write(target, "\tF");
        write(target, context.fixIdent[graphIdx]);
    });
    forEach(record.bppMatrGraphs, [&](RnaStructureGraph const & graph)
    {
        graphIdx = _findIndex(context.bppLabels, graph.specs);
        if (graphIdx == length(context.bppLabels))
            SEQAN_THROW(ParseError("ERROR: Bpp matrix description not found in header."));

        write(target, "\tM");
        write(target, context.bppIdent[graphIdx]);
    });
    writeValue(target, '\n');

    // write matrix
    unsigned offset = record.offset > 0 ? record.offset : 1;
    for (unsigned line = 0; line < record.seqLen; ++line)
    {
        // write index
        appendNumber(target, line + offset);

        // write sequence
        writeValue(target, '\t');
        write(target, record.sequence[line]);

        // write quality
        if (!empty(record.quality))
        {
            writeValue(target, '\t');
            writeValue(target, record.quality[line]);
        }

        // write reactivity
        forEach(record.typeID, [&](std::size_t tID)
        {
            if (tID >= length(record.reactivity) || tID >= length(record.reactError))
            {
                SEQAN_THROW(ParseError("ERROR: Invalid typeID " + std::to_string(tID) + " for printing reactivity."));
            }
            if (!empty(record.reactivity[tID]))
            {
                writeValue(target, '\t');
                appendNumber(target, record.reactivity[tID][line]);
            }
            if (!empty(record.reactError[tID]))
            {
                writeValue(target, '\t');
                appendNumber(target, record.reactError[tID][line]);
            }
        });

        // write fixed structure
        forEach(record.fixedGraphs, [&](RnaStructureGraph const & graph)
        {
            writeValue(target, '\t');
            if (degree(graph.inter, line) != 0)
            {
                RnaAdjacencyIterator adj_it(graph.inter, line);
                appendNumber(target, value(adj_it) + offset);
            }
            else
            {
                writeValue(target, '0');
            }
        });

        // write bpp matrix
        forEach(record.bppMatrGraphs, [&](RnaStructureGraph const & graph)
        {
            writeValue(target, '\t');
            if (degree(graph.inter, line) != 0)
            {
                RnaAdjacencyIterator adj_it(graph.inter, line);
                writeValue(target, '<');
                appendNumber(target, value(adj_it) + offset);
                writeValue(target, '/');
                appendNumber(target, cargo(findEdge(graph.inter, value(adj_it), line)));
                goNext(adj_it);
                while (!atEnd(adj_it))
                {
                    write(target, " | ");
                    appendNumber(target, value(adj_it) + offset);
                    writeValue(target, '/');
                    appendNumber(target, cargo(findEdge(graph.inter, value(adj_it), line)));
                    goNext(adj_it);
                }
                writeValue(target, '>');
            }
            else
            {
                write(target, "<>");
            }
        });
        writeValue(target, '\n');
    }
}

}  // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_EBPSEQ_READ_WRITE_H_
