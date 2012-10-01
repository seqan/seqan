// FRAGMENT(header)
#include <fstream>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/stream.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "ERROR: Invalid argument count." << std::endl
                  << "USAGE: " << argv[0] << " FILE" << std::endl;
        return 1;
    }
    
// FRAGMENT(read-sequences-single-pass)
    std::cerr << "Reading from std::fstream..." << std::endl << std::endl;
    // Open file and create RecordReader.
    std::fstream fasta(argv[1], std::ios_base::in | std::ios_base::binary);
    if (!fasta.good())
        return 1;
    RecordReader<std::fstream, SinglePass<> > reader(fasta);
    
    // Define variables for storing the sequences and sequence ids.
    StringSet<CharString> ids;
    StringSet<String<Dna5Q> > seqs;
    if (read2(ids, seqs, reader, Fastq()) != 0)
    {
        std::cerr << "ERROR reading FASTA." << std::endl;
        return 1;
    }
    
    // Write ids, sequences and qualities to output.
    typedef Iterator<StringSet<CharString>, Rooted>::Type TIdIter;
    typedef Iterator<StringSet<String<Dna5Q> >, Standard>::Type TSeqIter;
    TIdIter idIt = begin(ids, Rooted());
    TSeqIter seqIt = begin(seqs, Standard());
    CharString quals;
    for (; !atEnd(idIt); goNext(idIt), ++seqIt)
    {
        resize(quals, length(*seqIt));
        assignQualities(quals, *seqIt);
        std::cout << *idIt << '\t' << *seqIt << '\t' << quals << std::endl;
    }

// FRAGMENT(read-sequences-double-pass)
    std::cerr << std::endl << "Reading from memory mapped string..." << std::endl << std::endl;
    // Open file and create RecordReader.
    typedef String<char, MMap<> > TMMapString;
    TMMapString mmapString;
    if (!open(mmapString, argv[1], OPEN_RDONLY))
        return 1;
    RecordReader<TMMapString, DoublePass<Mapped> > reader2(mmapString);
    
    // Define variables for storing the sequences and sequence ids.
    StringSet<CharString> ids2;
    StringSet<String<Dna5Q> > seqs2;
    if (read2(ids2, seqs2, reader2, Fastq()) != 0)
    {
        std::cerr << "ERROR reading FASTA." << std::endl;
        return 1;
    }
    
    // Write ids, sequences and qualities to output.
    typedef Iterator<StringSet<CharString>, Rooted>::Type TIdIter;
    typedef Iterator<StringSet<String<Dna5Q> >, Standard>::Type TSeqIter;
    TIdIter idIt2 = begin(ids2, Rooted());
    TSeqIter seqIt2 = begin(seqs2, Standard());
    for (; !atEnd(idIt2); goNext(idIt2), ++seqIt2)
    {
        resize(quals, length(*seqIt2));
        assignQualities(quals, *seqIt2);
        std::cout << *idIt2 << '\t' << *seqIt2 << '\t' << quals << std::endl;
    }

// FRAGMENT(footer)    
    return 0;
}
