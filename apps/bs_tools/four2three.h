
#ifndef APPS_BS_TOOLS_FOUR2THREE_H_
#define APPS_BS_TOOLS_FOUR2THREE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fstream>

using namespace std;
using namespace seqan;


struct ConvertCT : public ::std::unary_function<char,char>
{
    inline char operator()(char x) const
    {
        if (('C' == x) || (x == 'c')) return ('T');
        return x;
    }
};

struct ConvertGA : public ::std::unary_function<char,char>
{
    inline char operator()(char x) const
    {
        if (('G' == x) || (x == 'g')) return ('A');
        return x;
    }
};

template<typename TOptions>
bool
preProcess(TOptions &options)
{
    typedef ModifiedString< CharString, ModView<ConvertCT> > TModCT;
    typedef ModifiedString< CharString, ModView<ConvertGA> > TModGA;

    CharString id, seq, quals;

    seqan::SeqFileIn  seqFileIn(toCString(options.inputFileName));
    seqan::SeqFileOut seqFileOut(toCString(options.outputFileName));

    // Get lower case of the output file name.  File endings are accepted in both upper and lower case.
    CharString tmp = options.inputFileName;
    toLower(tmp);

    if (endsWith(tmp, ".fa") || endsWith(tmp, ".fasta"))      // FASTA
    {
        while (!atEnd(seqFileIn))
        {
            readRecord(id, seq, seqFileIn);
            if (options.ctConversion)       // CT
            {
                TModCT modCT(seq);
                CharString tmp = modCT;
                writeRecord(seqFileOut, id, tmp);
            }
            else                            // GA
            {
                TModGA modGA(seq);
                CharString tmp = modGA;
                writeRecord(seqFileOut, id, tmp);
            }
        }
    }
    else if (endsWith(tmp, ".fastq") || endsWith(tmp, ".fq"))   // FASTQ
    {
        while (!atEnd(seqFileIn))
        {
            readRecord(id, seq, quals, seqFileIn);
            if (options.ctConversion)   //CT
            {
                TModCT modCT(seq);
                CharString tmp = modCT;
                writeRecord(seqFileOut, id, tmp, quals);
            }
            else                        // GA
            {
                TModGA modGA(seq);
                CharString tmp = modGA;
                writeRecord(seqFileOut, id, tmp, quals);
           }
        }
    }
    else
    {
       std::cerr << "ERROR: Something wrong with input!\n";
       return 1;
    }

    return 0;
}

#endif
