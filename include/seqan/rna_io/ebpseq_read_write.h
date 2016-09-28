// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Ebpseq
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Ebpseq
 * @headerfile <seqan/bpseq_io.h>
 * @brief Variant callinf format file.
 *
 * @signature typedef Tag<Ebpseq_> Ebpseq;
 */
struct Ebpseq_;
typedef Tag<Ebpseq_> Ebpseq;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [EbpseqHeader]
// ----------------------------------------------------------------------------

//inline void
//_parseEbpseqContig(CharString & chromName, CharString const & headerValue)
//{
//    if (length(headerValue) < 3u)
//        return;
//
//    CharString tmp = infix(headerValue, 1, length(headerValue) - 2);
//    StringSet<CharString> tmp2;
//    strSplit(tmp2, tmp, EqualsChar<','>());
//    for (unsigned i = 0; i < length(tmp2); ++i)
//    {
//        if (!startsWith(tmp2[i], "ID="))
//            continue;
//        chromName = suffix(tmp2[i], 3);
//    }
//}

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
        skipUntil(iter, IsWhitespace());
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(buffer, iter, IsNewline());

        switch (key)
        {
            case 'G': append(header.description, buffer); break;

            case 'S': appendValue(header.seqname, buffer); break;

            case 'F': appendValue(header.probmat, buffer); break;

            case 'M': appendValue(header.method, buffer); break;

            case 'T': appendValue(header.type, buffer); break;

            default: throw std::runtime_error("ERROR: Unknown key in header.");
        }

        clear(buffer);
        skipUntil(iter, NotFunctor<IsWhitespace>()); // move to next ##
    }
    header.amount = length(header.seqname);
    context.header = &header;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [EbpseqRecord]
// ----------------------------------------------------------------------------
// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter>
inline void
readRecord(RnaRecord & record, RnaIOContext & context, TForwardIter & iter, Ebpseq const & /*tag*/)
{
    std::string buffer;
    clear(record);

    // read sequence and type indices
    skipUntil(iter, NotFunctor<EqualsChar<'#'> >());
    while (value(iter) != '\n' && value(iter) != '\r')
    {
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(buffer, iter, IsWhitespace());
        if (startsWith(buffer, "S"))
        {
            if (lexicalCast(record.recordID, buffer.substr(1)))
                record.name = context.header->seqname[record.recordID - 1];
            else
                throw BadLexicalCast(record.recordID, buffer.substr(1));
        }
        else if (startsWith(buffer, "T"))
        {
            unsigned idx;
            if (lexicalCast(idx, buffer.substr(1)))
                appendValue(record.typeID, idx);
            else
                throw BadLexicalCast(idx, buffer.substr(1));
        }
        else
        {
            std::cerr << "Found tag '" << buffer << "'." << std::endl;
            throw std::runtime_error("ERROR: Unknown tag symbol.");
        }
        clear(buffer);
        skipUntil(iter, NotFunctor<IsBlank>());
    }
    if (record.recordID == record.UNDEF)
        throw std::runtime_error("ERROR: Record must contain a sequence tag.");

    // read column semantics
    StringSet<std::string> columnLabels {};
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

    // error checks and allocations
    //TODO check for Fi if i<len(fixedGraphs)
    unsigned line {1};
    resize(record.reactivity, length(record.typeID));
    resize(record.reactError, length(record.typeID));
    resize(record.fixedGraphs, length(context.header->method));
    resize(record.bppMatrGraphs, length(context.header->probmat));

    // read matrix
    while (!atEnd(iter) && value(iter) != '#')
    {
        for (std::string const label : columnLabels)
        {
            unsigned idx, pairPos;
            readUntil(buffer, iter, IsWhitespace());

            if (label == "I")
            {
                if (line == 1u && !lexicalCast(record.offset, buffer))
                    throw BadLexicalCast(record.offset, buffer);
            }
            else if (label == "NT") {
                appendValue(record.sequence, buffer[0]);
            }
            else if (label == "QU") {
                appendValue(record.quality, buffer[0]);
            }
            else if (startsWith(label, "F") || startsWith(label, "M")) {
                bool fixed = startsWith(label, "F");
                lexicalCast(idx, label.substr(1));
                TRnaRecordGraph & graph = fixed ? record.fixedGraphs[idx - 1].inter
                                                : record.bppMatrGraphs[idx - 1].inter;
                addVertex(graph);
                if (fixed)
                {
                    if (!lexicalCast(pairPos, buffer))
                        throw BadLexicalCast(pairPos, buffer);
                    if (pairPos != 0 && line + record.offset - 1 > pairPos)    // add edge if base is connected
                        addEdge(graph, pairPos - record.offset, line - 1, 1.);
                }
                else
                {
                    if (buffer.find('<') == std::string::npos)
                        throw std::runtime_error("ERROR: Expected a bracket pair: < >");
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
                            throw BadLexicalCast(pairPos, match.str().substr(0, pos));
                        if (!lexicalCast(prob, match.str().substr(pos + 1)))
                            throw BadLexicalCast(prob, match.str().substr(pos + 1));
                        if (pairPos != 0 && line + record.offset - 1 > pairPos)
                            addEdge(graph, pairPos - record.offset, line - 1, prob);
                        buffer = match.suffix();
                    }
                }
            }
            else if (startsWith(label, "RE"))
            {
                lexicalCast(idx, label.substr(2));
                unsigned typeID {0};
                while (record.typeID[typeID] != idx)
                    ++typeID;
                float react;
                if (!lexicalCast(react, buffer))
                    throw BadLexicalCast(react, buffer);
                appendValue(record.reactError[typeID], react);
            }
            else if (startsWith(label, "R"))
            {
                lexicalCast(idx, label.substr(1));
                unsigned typeID {0};
                while (record.typeID[typeID] != idx)
                    ++typeID;
                float react;
                if (!lexicalCast(react, buffer))
                    throw BadLexicalCast(react, buffer);
                appendValue(record.reactivity[typeID], react);
            }
            else
            {
                std::cerr << "WARNING: label '" << label << "' not implemented." << std::endl;
            }

            skipUntil(iter, NotFunctor<IsBlank>());
            clear(buffer);
        }
        ++line;
        skipLine(iter);
    }
    record.seqLen = length(record.sequence);
    appendValue(context.record, &record);
}

// ----------------------------------------------------------------------------
// Function writeHeader()                                           [EbpseqHeader]
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

    for (unsigned idx = 1; idx <= length(header.seqname); ++idx)
    {
        write(target, "## S");
        write(target, idx);
        write(target, ": ");
        write(target, header.seqname[idx-1]);
        writeValue(target, '\n');
    }

    for (unsigned idx = 1; idx <= length(header.method); ++idx)
    {
        write(target, "## F");
        write(target, idx);
        write(target, ": ");
        write(target, header.method[idx-1]);
        writeValue(target, '\n');
    }

    for (unsigned idx = 1; idx <= length(header.probmat); ++idx)
    {
        write(target, "## M");
        write(target, idx);
        write(target, ": ");
        write(target, header.probmat[idx-1]);
        writeValue(target, '\n');
    }

    for (unsigned idx = 1; idx <= length(header.type); ++idx)
    {
        write(target, "## T");
        write(target, idx);
        write(target, ": ");
        write(target, header.type[idx-1]);
        writeValue(target, '\n');
    }
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [EbpseqRecord]
// ----------------------------------------------------------------------------

template <typename TTarget>
inline void
writeRecord(TTarget & target, RnaRecord const & record, RnaIOContext & context, Ebpseq const & /*tag*/)
{
    unsigned idx;
    // write head line with S and T tags
    write(target, "# S");
    write(target, record.recordID);
    for (idx = 0; idx < length(record.typeID); ++idx)
    {
        write(target, " T");
        write(target, record.typeID[idx]);
    }
    writeValue(target, '\n');

    // write column header
    write(target, "# I\tNT");
    if (!empty(record.quality))
        write(target, "\tQU");
    for (idx = 0; idx < length(record.typeID); ++idx)
    {
        if (!empty(record.reactivity[idx]))
        {
            write(target, "\tR");
            write(target, record.typeID[idx]);
        }
        if (!empty(record.reactError[idx]))
        {
            write(target, "\tRE");
            write(target, record.typeID[idx]);
        }
    }
    for (idx = 0; idx < length(record.fixedGraphs); ++idx)
    {
        if (numVertices(record.fixedGraphs[idx].inter) > 0)
        {
            write(target, "\tF");
            write(target, idx + 1);
        }
    }
    for (idx = 0; idx < length(record.bppMatrGraphs); ++idx)
    {
        if (numVertices(record.bppMatrGraphs[idx].inter) > 0)
        {
            write(target, "\tM");
            write(target, idx + 1);
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
            if (!empty(record.reactivity[idx]))
            {
                writeValue(target, '\t');
                write(target, record.reactivity[idx][line]);
            }
            if (!empty(record.reactError[idx]))
            {
                writeValue(target, '\t');
                write(target, record.reactError[idx][line]);
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
                    TRnaAdjacencyIterator adj_it(record.fixedGraphs[idx].inter, line);
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
            TRnaRecordGraph const & graph = record.bppMatrGraphs[idx].inter;
            if (numVertices(graph) > 0)
            {
                writeValue(target, '\t');
                if (degree(graph, line) != 0)
                {
                    TRnaAdjacencyIterator adj_it(graph, line);
                    writeValue(target, '<');
                    write(target, value(adj_it) + offset);
                    writeValue(target, '/');
                    write(target, cargo(findEdge(graph, value(adj_it), line)));
                    goNext(adj_it);
                    while (!atEnd(adj_it))
                    {
                        write(target, " | ");
                        write(target, value(adj_it) + offset);
                        writeValue(target, '/');
                        write(target, cargo(findEdge(graph, value(adj_it), line)));
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
