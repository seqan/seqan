// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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

inline unsigned _findIndex(StringSet<CharString> const & set, CharString const & str)
{
    unsigned idx = 0;
    while (idx < length(set) && set[idx] != str)
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
    std::string buffer;
    unsigned idx;
    clear(record);

    // read sequence and type indices
    skipUntil(iter, NotFunctor<EqualsChar<'#'> >());
    while (value(iter) != '\n' && value(iter) != '\r')
    {
        bool hadErr = false;
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(buffer, iter, IsWhitespace());
        if (startsWith(buffer, "S"))
        {
            idx = _findIndex(context.seqIdent, buffer.substr(1));
            record.recordID = idx;
            if (idx < length(context.seqIdent))
            {
                record.name = context.seqLabels[idx];
            }
            else
            {
                hadErr = true;
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
                hadErr = true;
            }
        }
        else
        {
            hadErr = true;
        }
        if (hadErr)
        {
            std::cerr << "Found tag '" << buffer << "'." << std::endl;
            SEQAN_THROW(ParseError("ERROR: Unknown tag symbol."));
        }
        clear(buffer);
        skipUntil(iter, NotFunctor<IsBlank>());
    }
    if (record.hasUndefinedID())
    {
        SEQAN_THROW(ParseError("ERROR: Record must contain a sequence tag."));
    }

    // read column semantics
    StringSet<std::string> columnLabels{};
    skipUntil(iter, NotFunctor<IsWhitespace>());
    readOne(buffer, iter, EqualsChar<'#'>());
    clear(buffer);
    while (value(iter) != '\n' && value(iter) != '\r')
    {
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(buffer, iter, IsWhitespace());
        appendValue(columnLabels, buffer);
        clear(buffer);
        skipUntil(iter, NotFunctor<IsBlank>());
    }
    skipUntil(iter, NotFunctor<IsWhitespace>());

    // allocations
    unsigned line {1};
    resize(record.reactivity, length(context.typIdent));
    resize(record.reactError, length(context.typIdent));
    resize(record.fixedGraphs, length(context.fixIdent));
    resize(record.bppMatrGraphs, length(context.bppIdent));

    // read matrix
    while (!atEnd(iter) && value(iter) != '#')
    {
        for (std::string const label : columnLabels)
        {
            bool hadErr = false;
            unsigned pairPos;
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
                idx = _findIndex(context.fixIdent, label.substr(1));
                if (idx < length(context.fixIdent))
                {
                    RnaStructureGraph & graph = record.fixedGraphs[idx];
                    addVertex(graph.inter);
                    graph.specs = context.fixLabels[idx];
                    if (!lexicalCast(pairPos, buffer))
                    {
                        SEQAN_THROW(BadLexicalCast(pairPos, buffer));
                    }
                    if (pairPos != 0 && line + record.offset - 1 > pairPos)    // add edge if base is connected
                    {
                        addEdge(graph.inter, pairPos - record.offset, line - 1, 1.0);
                    }
                }
                else
                {
                    hadErr = true;
                }
            }
            else if (startsWith(label, "M"))
            {
                idx = _findIndex(context.bppIdent, label.substr(1));
                if (idx < length(context.bppIdent))
                {
                    RnaStructureGraph & graph = record.bppMatrGraphs[idx];
                    addVertex(graph.inter);
                    graph.specs = context.bppLabels[idx];
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
                else
                {
                    hadErr = true;
                }
            }
            else if (startsWith(label, "RE"))
            {
                idx = _findIndex(context.typIdent, label.substr(2));
                if (idx < length(context.typIdent))
                {
                    float react;
                    if (!lexicalCast(react, buffer))
                    {
                        SEQAN_THROW(BadLexicalCast(react, buffer));
                    }
                    appendValue(record.reactError[idx], react);
                }
                else
                {
                    hadErr = true;
                }
            }
            else if (startsWith(label, "R"))
            {
                idx = _findIndex(context.typIdent, label.substr(1));
                if (idx < length(context.typIdent))
                {
                    float react;
                    if (!lexicalCast(react, buffer))
                    {
                        SEQAN_THROW(BadLexicalCast(react, buffer));
                    }
                    appendValue(record.reactivity[idx], react);
                }
                else
                {
                    hadErr = true;
                }
            }
            else
            {
                hadErr = true;
            }
            if (hadErr)
            {
                std::cerr << "Found label '" << label << "'." << std::endl;
                SEQAN_THROW(ParseError("ERROR: Unknown label."));
            }
            skipUntil(iter, NotFunctor<IsBlank>());
            clear(buffer);
        }
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
    clear(context);
    if (!empty(header.description))
    {
        write(target, "## G: ");
        write(target, header.description);
        writeValue(target, '\n');
    }

    for (unsigned idx = 1; idx <= length(header.seqLabels); ++idx)
    {
        write(target, "## S");
        write(target, idx);
        write(target, ": ");
        write(target, header.seqLabels[idx - 1]);
        writeValue(target, '\n');
        appendValue(context.seqIdent, std::to_string(idx));
    }

    for (unsigned idx = 1; idx <= length(header.fixLabels); ++idx)
    {
        write(target, "## F");
        write(target, idx);
        write(target, ": ");
        write(target, header.fixLabels[idx - 1]);
        writeValue(target, '\n');
        appendValue(context.fixIdent, std::to_string(idx));
    }

    for (unsigned idx = 1; idx <= length(header.bppLabels); ++idx)
    {
        write(target, "## M");
        write(target, idx);
        write(target, ": ");
        write(target, header.bppLabels[idx - 1]);
        writeValue(target, '\n');
        appendValue(context.bppIdent, std::to_string(idx));
    }

    for (unsigned idx = 1; idx <= length(header.typeLabels); ++idx)
    {
        write(target, "## T");
        write(target, idx);
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
    std::size_t numFix{}, numBpp{};
    for (unsigned idx = 0; idx < length(records); ++idx)
    {
        records[idx].recordID = idx;
        if (!empty(records[idx].comment))
        {
            if (!empty(header.description))
            {
                append(header.description, "\t");
            }
            append(header.description, records[idx].comment);
        }

        appendValue(header.seqLabels, records[idx].name);

        if (length(records[idx].fixedGraphs) > numFix)
        {
            numFix = length(records[idx].fixedGraphs);
        }
        if (length(records[idx].bppMatrGraphs) > numBpp)
        {
            numBpp = length(records[idx].bppMatrGraphs);
        }
    }

    for (unsigned idx = 1; idx <= numFix; ++idx)
        appendValue(header.fixLabels, "n/a");

    for (unsigned idx = 1; idx <= numBpp; ++idx)
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

    unsigned idx;
    // write head line with S and T tags
    write(target, "# S");
    write(target, context.seqIdent[record.recordID]);
    for (idx = 0; idx < length(record.typeID); ++idx)
    {
        if (record.typeID[idx] >= length(context.typIdent))
        {
            SEQAN_THROW(ParseError("ERROR: Type ID exceeds number of types in ebpseq header."));
        }
        write(target, " T");
        write(target, context.typIdent[record.typeID[idx]]);
    }
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
    for (idx = 0; idx < length(record.typeID); ++idx)
    {
        if (!empty(record.reactivity[record.typeID[idx]]))
        {
            write(target, "\tR");
            write(target, context.typIdent[record.typeID[idx]]);
        }
        if (!empty(record.reactError[record.typeID[idx]]))
        {
            write(target, "\tRE");
            write(target, context.typIdent[record.typeID[idx]]);
        }
    }
    for (idx = 0; idx < length(record.fixedGraphs); ++idx)
    {
        if (numVertices(record.fixedGraphs[idx].inter) > 0)
        {
            write(target, "\tF");
            write(target, context.fixIdent[idx]);
        }
    }
    for (idx = 0; idx < length(record.bppMatrGraphs); ++idx)
    {
        if (numVertices(record.bppMatrGraphs[idx].inter) > 0)
        {
            write(target, "\tM");
            write(target, context.bppIdent[idx]);
        }
    }
    writeValue(target, '\n');

    // write matrix
    unsigned offset = record.offset > 0 ? record.offset : 1;
    for (unsigned line = 0; line < record.seqLen; ++line)
    {
        // write index
        write(target, line + offset);

        // write sequence
        writeValue(target, '\t');
        write(target, record.sequence[line]);

        // write quality
        if (!empty(record.quality))
        {
            writeValue(target, '\t');
            write(target, record.quality[line]);
        }

        // write reactivity
        for (idx = 0; idx < length(record.typeID); ++idx)
        {
            if (record.typeID[idx] >= length(record.reactivity) || record.typeID[idx] >= length(record.reactError))
            {
                SEQAN_THROW(ParseError("ERROR: Invalid typeID for printing reactivity."));
            }
            if (!empty(record.reactivity[record.typeID[idx]]))
            {
                writeValue(target, '\t');
                write(target, record.reactivity[record.typeID[idx]][line]);
            }
            if (!empty(record.reactError[record.typeID[idx]]))
            {
                writeValue(target, '\t');
                write(target, record.reactError[record.typeID[idx]][line]);
            }
        }

        // write fixed structure
        for (idx = 0; idx < length(record.fixedGraphs); ++idx)
        {
            if (numVertices(record.fixedGraphs[idx].inter) > 0)
            {
                writeValue(target, '\t');
                if (degree(record.fixedGraphs[idx].inter, line) != 0)
                {
                    RnaAdjacencyIterator adj_it(record.fixedGraphs[idx].inter, line);
                    write(target, value(adj_it) + offset);
                }
                else
                {
                    writeValue(target, '0');
                }
            }
        }

        // write bpp matrix
        for (idx = 0; idx < length(record.bppMatrGraphs); ++idx)
        {
            RnaStructureGraph const & graph = record.bppMatrGraphs[idx];
            if (numVertices(graph.inter) > 0)
            {
                writeValue(target, '\t');
                if (degree(graph.inter, line) != 0)
                {
                    RnaAdjacencyIterator adj_it(graph.inter, line);
                    writeValue(target, '<');
                    write(target, value(adj_it) + offset);
                    writeValue(target, '/');
                    write(target, cargo(findEdge(graph.inter, value(adj_it), line)));
                    goNext(adj_it);
                    while (!atEnd(adj_it))
                    {
                        write(target, " | ");
                        write(target, value(adj_it) + offset);
                        writeValue(target, '/');
                        write(target, cargo(findEdge(graph.inter, value(adj_it), line)));
                        goNext(adj_it);
                    }
                    writeValue(target, '>');
                }
                else
                {
                    write(target, "<>");
                }
            }
        }
        writeValue(target, '\n');
    }
}

}  // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_RNA_IO_EBPSEQ_READ_WRITE_H_
