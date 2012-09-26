///A tutorial about the use of files and file formats.
#include <iostream>
#include <fstream>
#include <cstdio>

#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
///Open a standard library stream for binary write. 
///You can use both C++ iostreams or old style FILE pointer.
///The file format is Fasta.
	std::FILE * fl = std::fopen("testfile.fa", "wb");
    write(fl, "aacagtattagaccactaggaccct", "a test file", Fasta());
	close (fl);
///Read a fasta file.
	std::fstream fstrm;
	fstrm.open("testfile.fa", std::ios_base::in | std::ios_base::binary);
	String<char> fasta_tag;
	String<Dna> fasta_seq;
///Read the meta-information.
	readMeta(fstrm, fasta_tag, Fasta());
	std::cout << fasta_tag << "\n";	//prints "a test file"
///Read the sequence.
	read(fstrm, fasta_seq, Fasta());
	std::cout << fasta_seq << "\n";	//prints the sequence
	fstrm.close();
///Open a file using a file reader string.
	String<Dna, FileReader<Fasta> > fr("testfile.fa");
	std::cout << fr << "\n";			//prints the sequence
	return 0;
}

