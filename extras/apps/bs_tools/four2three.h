
#ifndef CORE_APPS_BS_TOOLS_FOUR2THREE_H_
#define CORE_APPS_BS_TOOLS_FOUR2THREE_H_

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

    CharString id;
    CharString seq;
    CharString quals;
   

    SequenceStream streamRead(toCString(options.inputFileName));
    SequenceStream streamWrite(toCString(options.outputFileName), SequenceStream::WRITE);
   
    // Get lower case of the output file name.  File endings are accepted in both upper and lower case.
    CharString tmp = options.inputFileName;
    toLower(tmp);

    if (endsWith(tmp, ".fa") || endsWith(tmp, ".fasta"))      // FASTA
    {
        if (!isGood(streamRead))
        {
            std::cerr << "ERROR: Could not open ref file.\n";
            return 1;
        }
        while (!atEnd(streamRead))
        {
            if (readRecord(id, seq, streamRead) != 0)                                
            {
                std::cerr << "ERROR: Could not read from reference file!\n";
                return 1;
            }
            if (options.ctConversion)       // CT
            {
                TModCT modCT(seq);
                CharString tmp = modCT;
                if (writeRecord(streamWrite, id, tmp) != 0)
                {
                    std::cerr << "ERROR: Could not write to file!\n";
                    return 1;
                }
            }
            else                            // GA
            {
                TModGA modGA(seq);
                CharString tmp = modGA;                                             
                if (writeRecord(streamWrite, id, tmp) != 0)
                {
                    std::cerr << "ERROR: Could not write to file!\n";
                    return 1;
                }
            }
        }
    }
    else if (endsWith(tmp, ".fastq") || endsWith(tmp, ".fq"))   // FASTQ
    {
        if (!isGood(streamRead))
        {
            std::cerr << "ERROR: Could not open read file.\n";
            return 1;
        }
        while (!atEnd(streamRead))
        {
            if (readRecord(id, seq, quals, streamRead) != 0)
            {
                std::cerr << "ERROR: Could not read from read file!\n";
                return 1;
            }
            if (options.ctConversion)   //CT
            {
                TModCT modCT(seq);
                CharString tmp = modCT;
                if (writeRecord(streamWrite, id, tmp, quals) != 0)
                {
                    std::cerr << "ERROR: Could not write to file!\n";
                    return 1;
                }
            }
            else                        // GA
            { 
                TModGA modGA(seq);
                CharString tmp = modGA; 
                 if (writeRecord(streamWrite, id, tmp, quals) != 0)
                {
                    std::cerr << "ERROR: Could not write to file!\n";
                    return 1;
                }
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
