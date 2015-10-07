///A tutorial about the use of files and file formats.
#include <iostream>
#include <fstream>
#include <cstdio>

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
///Open a standard library stream for binary write.
///You can use both C++ iostreams or old style FILE pointer.
///The file format is Fasta.
    std::ofstream fl("testfile.fa");
    DirectionIterator<std::ofstream, Output>::Type out = directionIterator(fl, Output());
    writeRecord(out, "a test file", "aacagtattagaccactaggaccct", Fasta());
    fl.close();
///Read a fasta file.
    std::ifstream fstrm("testfile.fa");
    DirectionIterator<std::ifstream, Input>::Type in = directionIterator(fstrm, Input());
    String<char> fasta_tag;
    String<Dna> fasta_seq;
///Read the meta-information and the sequence.
    readRecord(fasta_tag, fasta_seq, in, Fasta());
    std::cout << fasta_tag << "\n"; //prints "a test file"
    std::cout << fasta_seq << "\n"; //prints the sequence
    fstrm.close();
///Open a file using a file reader string.
    SeqFileIn sfin("testfile.fa");
    readRecord(fasta_tag, fasta_seq, sfin);
    std::cout << fasta_seq << "\n";         //prints the sequence
    return 0;
}
