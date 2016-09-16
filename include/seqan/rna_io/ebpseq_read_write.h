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
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_READ_WRITE_BPSEQ_H_
#define SEQAN_INCLUDE_SEQAN_BPSEQ_READ_WRITE_BPSEQ_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Bpseq
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Bpseq
 * @headerfile <seqan/bpseq_io.h>
 * @brief Variant callinf format file.
 *
 * @signature typedef Tag<Bpseq_> Bpseq;
 */
struct Bpseq_;
typedef Tag<Bpseq_> Bpseq;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BpseqHeader]
// ----------------------------------------------------------------------------

//inline void
//_parseBpseqContig(CharString & chromName, CharString const & headerValue)
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

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readHeader(RnaHeader & header,
           BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Bpseq const & /*tag*/)
{
    clear(header);
    CharString buffer;
    RnaHeaderRecord record;


    clear(buffer);
    skipOne(iter);
    if(value(iter) == '#')
    {
        // Is header line.
        skipUntil(iter, NotFunctor<IsWhitespace>());
        clear(record);

        // Read header key.
        readUntil(record.key, iter, EqualsChar<':'>());

        // Parse out names of sequences
        while (!startsWith(record.key, "S")
        {
            // Skip ':'.  
            skipOne(iter);

            // Read header value.
            readLine(record.value, iter);
            appendValue(header, record);
            
            appendName(contigNamesCache(context), record.value);
            skipOne(iter); //skip to newline
            skipOne(iter); //skip first '#'
            if(value(iter)!= '#')
                break;  ///SHOULD THIS BE RETURN BECAUSE THEN WE ARE DONE WITH HEADER INFO IF THERE ARE NO MORE ## LINES?
            else
                readUntil(record.key, iter, EqualsChar<':'>());
        }

        while(startsWith(record.key, "F")){
            // Skip ':'.  
            skipOne(iter);

            // Read header value.
            readLine(record.value, iter);
            appendValue(header, record);
            
            appendName(fixedStructureNamesCache(context), record.value);
            skipOne(iter); //skip to newline

            //Read new fixed structure key
            if(value(iter)!= '#')
                break;  
            else
                readUntil(record.key, iter, EqualsChar<':'>());
        }

        while(startsWith(record.key, "M"))
        {
            // Skip ':'.  
            skipOne(iter);

            // Read header value.
            readLine(record.value, iter);
            appendValue(header, record);
            
            appendName(baseProbabilityNamesCache(context), record.value);
            skipOne(iter); //skip to newline

            //Read new base probability key
            if(value(iter)!= '#')
                break;  
            else
                readUntil(record.key, iter, EqualsChar<':'>());
        }
    }


    appendValue(header, record);

    
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BpseqRecord]
// ----------------------------------------------------------------------------
// Read record, updating list of known sequences if new one occurs.

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(RnaRecord & record,
           BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Bpseq const & /*tag*/)
{
    typedef OrFunctor<IsSpace, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bpseq> > NextEntry;
    clear(record);
    CharString &buffer = context.buffer;
    CharString tmpStr="";
    unsigned counter = 0;

    //temporary holders until I figure out a better way to do this
    int quality_flag;
    int reaction_1_flag;
    int reaction_2_flag;
    int reaction_error_1_flag;
    int reaction_error_2_flag;
    String<int> fixed_structure_flags;
    //resize to be size specified in header/context
    String<int> matrix_probablilty_pair_flag;
    //resize to be size specificed in header/context

    //FOR T1 and T2, save in an array and set corresponding info to corresponding records

    if(value(iter)== '#'){
        skipOne(iter);
        while(value(iter) != '\n'){
            skipOne(iter);
            readUntil(context.buffer, iter, IsWhitespace());
            clear(context.buffer);
            //Check with context, compare and store info based on name or whatever
            //NEED A FUNCTION TO DO THAT
        }
        skipOne(iter);//Skip to new line, header of the info in the columns
        skipOne(iter); //skip #
        skipOne(iter);  //skip until start
        readUntil(context.buffer, iter, IsWhitespace());

        //INDEX
        if(context.buffer != "I"){
            std::cerr << "ERROR: Incorrect column information" << std::endl;   
            return;
        }
        skipOne(iter);
        clear(context.buffer);

        //NT
        readUntil(context.buffer, iter, IsWhitespace());
        if(context.buffer  != "NT"){
            std::cerr << "ERROR: Incorrect column information" << std::endl;   
            return;
        }
        skipOne(iter);
        clear(context.buffer);

        //QUALITY
        readUntil(context.buffer, iter, IsWhitespace());
        if(context.buffer  != "QU"){
            quality_flag = -1;
        }
        else{
            quality_flag = 1;
            skipOne(iter);
            clear(context.buffer);
            readUntil(context.buffer, iter, IsWhitespace());
        }

        //REACTION____1
        if(context.buffer != "R1"){
            reaction_1_flag = -1;
        }
        else{
            reaction_1_flag = 1;
            skipOne(iter);
            clear(context.buffer);
            readUntil(context.buffer, iter, IsWhitespace());
        }

        //REACTION 2
        if(context.buffer  != "R2"){
            reaction_2_flag = -1;
        }
        else{
            reaction_2_flag = 1;
            skipOne(iter);
            clear(context.buffer);
            readUntil(context.buffer, iter, IsWhitespace());
        }

        //REACTION ERROR 1
        if(context.buffer  != "RE1"){
            reaction_error_1_flag = -1;
        }
        else{
            reaction_error_1_flag = 1;
            skipOne(iter);
            clear(context.buffer);
            readUntil(context.buffer, iter, IsWhitespace());
        }

        //REACTION ERROR 2
        if(context.buffer  != "RE2"){
            reaction_error_2_flag = -1;
        }
        else{
            reaction_error_2_flag = 1;
            skipOne(iter);
            clear(context.buffer);
            readUntil(context.buffer, iter, IsWhitespace());
        }

        //A for loop to iterate through all F1, F2...FN sections
        String f_name;
        int f_pos;
        for(unsigned i = 0; i < /*Length specified in context of F1*/; ++i)
        {
            f_pos = i+1;
            f_name = "F" + f_pos;
            /////////TO-DO: check name and that this is correct way of doing this
            if(context.buffer != f_name){
                fixed_structure_flags[i] = -1;
            }
            else{
                fixed_structure_flags[i] = 1;
                skipOne(iter);
                clear(context.buffer);
                readUntil(context.buffer, iter, IsWhitespace());
            }
        }

        //A for loop to iterate through all M1, M2...MN sections
        String m_name;
        int m_pos;
        for(unsigned i = 0; i < /*Length specified in context of F1*/; ++i)
        {
            m_pos = i+1;
            m_name = "M" + f_pos;
            /////////TO-DO: check name and that this is correct way of doing this
            if(context.buffer != f_name){
                matrix_probablilty_pair_flag[i] = -1;
            }
            else{
                matrix_probablilty_pair_flag[i] = 1;
                skipOne(iter);
                clear(context.buffer);
                readUntil(context.buffer, iter, IsWhitespace());
            }
        }

    }
    else
        std::cerr << "PARSE ERROR\n";
    //FINISHED WITH HEADER INFO

    while (!atEnd(iter) && value(iter) != '#')
    {
        // SEQPOS
        clear(buffer);
        readUntil(buffer, iter, NextEntry());
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("SEQPOS"));
//        record.seqPos = nameToId(contigNamesCache(context), buffer);
//        std::cout << lexicalCast<__int32>(buffer) << std::endl;//<< "buffer size = " << length(buffer);
//        ciao = lexicalCast<__int32>(buffer);
        if(counter == 0)
            record.begPos = lexicalCast<__int32>(buffer);
        skipUntil(iter, NextEntry());

        // SEQ
        clear(buffer);
        readUntil(buffer, iter, NextEntry());
        appendValue(record.sequence,buffer[0]);
        if (empty(record.seq))
            SEQAN_THROW(EmptyFieldError("SEQ"));

        skipUntil(iter, NextEntry());

        // PAIR
        clear(buffer);
        readUntil(buffer, iter, OrFunctor<IsSpace, IsNewline>());
        if (empty(buffer))
            SEQAN_THROW(EmptyFieldError("PAIR"));
//        record.pair = nameToId(contigNamesCache(context), buffer);
        appendValue(record.pair, lexicalCast<unsigned>(buffer));
//        skipOne(iter);
//        // the following columns are optional
//        if (atEnd(iter) || IsNewline()(value(iter)))
//        {
//            skipLine(iter);
//        }
        // skip line break and optional additional columns
        skipLine(iter);
        counter++;
    }
    if(record.begPos != 1)      //set beginning record position
        record.endPos = counter - record.begPos + 1;
    else
        record.endPos = counter;    //set end record position
    record.amount = record.endPos - record.begPos + 1;  //set amount of records

//    std::cout << "value(iter)" << value(iter) << std::endl;

    // CHROM
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("CHROM"));
    record.rID = nameToId(contigNamesCache(context), buffer);
    skipOne(iter);

    // POS
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("POS"));
    record.beginPos = lexicalCast<__int32>(buffer) - 1; // Translate from 1-based to 0-based.
    skipOne(iter);

    // ID
    readUntil(record.id, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("ID"));
    skipOne(iter);

    // REF
    readUntil(record.ref, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("REF"));
    skipOne(iter);

    // ALT
    readUntil(record.alt, iter, NextEntry());
    if (empty(record.id))
        SEQAN_THROW(EmptyFieldError("ALT"));
    skipOne(iter);

    // QUAL
    clear(buffer);
    readUntil(buffer, iter, NextEntry());
    if (empty(buffer))
        SEQAN_THROW(EmptyFieldError("QUAL"));

    if (buffer == ".")
        record.qual = BpseqRecord::MISSING_QUAL();
    else
        lexicalCastWithException(record.qual, buffer);

    skipOne(iter);

    // FILTER
    readUntil(record.filter, iter, NextEntry());
    if (empty(record.filter))
        SEQAN_THROW(EmptyFieldError("FILTER"));
    skipOne(iter);

    // INFO
    readUntil(record.info, iter, OrFunctor<IsTab, IsNewline>());
    if (empty(record.info))
        SEQAN_THROW(EmptyFieldError("INFO"));

    // the following columns are optional
    if (atEnd(iter) || IsNewline()(value(iter)))
    {
        skipLine(iter);
        return;
    }
    skipOne(iter);

    // FORMAT
    readUntil(record.format, iter, NextEntry());
    if (empty(record.format))
        SEQAN_THROW(EmptyFieldError("FORMAT"));
    skipOne(iter);

    // The samples.
    unsigned numSamples = length(sampleNames(context));
    for (unsigned i = 0; i < numSamples; ++i)
    {
        clear(buffer);
        if (i + 1 != numSamples)
        {
            readUntil(buffer, iter, NextEntry());
            skipOne(iter);
        }
        else
        {
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
        }

        if (empty(buffer))
        {
            char buffer[30];    // == 9 (GENOTYPE_) + 20 (#digits in MIN_INT64) + 1 (trailing zero)
            sprintf(buffer, "GENOTYPE_%u", i + 1);
            SEQAN_THROW(EmptyFieldError(buffer));
        }
        appendValue(record.genotypeInfos, buffer);
    }

    // skip line break and optional additional columns
    skipLine(iter);
*/
    return;
}

// ----------------------------------------------------------------------------
// Function writeHeader()                                           [BpseqHeader]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
writeHeader(TTarget & target,
            BpseqHeader const & header,
            BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            Bpseq const & /*tag*/)
{
    for (unsigned i = 0; i < length(header); ++i)
    {
        write(target, "## ");
        write(target, header[i].key);
        write(target, ': ');
        write(target, header[i].value);
        writeValue(target, '\n');
    }
//
//    write(target, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
//    for (unsigned i = 0; i < length(sampleNames(context)); ++i)
//    {
//        writeValue(target, '\t');
//        write(target, sampleNames(context)[i]);
//    }
//    writeValue(target, '\n');
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [BpseqRecord]
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
writeRecord(TTarget & target,
            BpseqRecord const & record,
            BpseqIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            Bpseq const & /*tag*/)
{
    write(target, "# ");
    write(target, "I\tNT\t");
    //THE FLAGS NEED TO BE STORED IN CONTEXT BUFFER!!
    if(quality_flag == 1)
        write(target, "QU\t")
    if(reaction_1_flag == 1)
        write(target, "R1\t")
    if(reaction_2_flag == 1)
        write(target, "R2\t")
    if(reaction_error_1_flag == 1)
        write(target, "RE1\t")
    if(reaction_error_2_flag == 1)
        write(target, "RE2\t")
    //F1...Fn and M1...Mn
    write(target, '\n');

    
//  // SEQPOS
//    for (unsigned i = 0; i < length(record.seqPos); ++i)
//    {
//        writeValue(target, '\t');
////        if (empty(record.genotypeInfos[i]))
////            writeValue(target, '.');
////        else
//            write(target, record.seqPos[i]);
//    }
////    write(target, record.seqPos);
//    writeValue(target, '\n');
//
//    // SEQ
//    write(target, record.seq);
//    writeValue(target, '\n');
//
//    // INTERPOS
//    for (unsigned i = 0; i < length(record.interPos); ++i)
//    {
//        writeValue(target, '\t');
////        if (empty(record.genotypeInfos[i]))
////            writeValue(target, '.');
////        else
//            write(target, record.interPos[i]);
//    }
////    write(target, record.interPos);
//    writeValue(target, '\n');

    // SEQPOS
    for (unsigned i = 0; i < length(record.seqPos); ++i)
    {
        write(target, record.seqPos[i]);
        writeValue(target, ' ');
        write(target, record.seq[i]);
        writeValue(target, ' ');
        write(target, record.interPos[i]);
        writeValue(target, '\n');
    }
/*
    // CHROM
    write(target, contigNames(context)[record.rID]);
    writeValue(target, '\t');

    // POS
    appendNumber(target, record.beginPos + 1);
    writeValue(target, '\t');

    // ID
    if (empty(record.id))
        writeValue(target, '.');
    else
        write(target, record.id);
    writeValue(target, '\t');

    // REF
    if (empty(record.ref))
        writeValue(target, '.');
    else
        write(target, record.ref);
    writeValue(target, '\t');

    // ALT
    if (empty(record.alt))
        writeValue(target, '.');
    else
        write(target, record.alt);
    writeValue(target, '\t');

    // QUAL
    if (record.qual != record.qual)  // only way to test for nan
        writeValue(target, '.');
    else
        appendNumber(target, record.qual);

    // FILTER
    writeValue(target, '\t');
    if (empty(record.filter))
        writeValue(target, '.');
    else
        write(target, record.filter);
    writeValue(target, '\t');

    // INFO
    if (empty(record.info))
        writeValue(target, '.');
    else
        write(target, record.info);

    // FORMAT
    writeValue(target, '\t');
    if (empty(record.format))
        writeValue(target, '.');
    else
        write(target, record.format);

    // The samples.
    for (unsigned i = 0; i < length(record.genotypeInfos); ++i)
    {
        writeValue(target, '\t');
        if (empty(record.genotypeInfos[i]))
            writeValue(target, '.');
        else
            write(target, record.genotypeInfos[i]);
    }
    writeValue(target, '\n');
*/
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BPSEQ_READ_WRITE_BPSEQ_H_
