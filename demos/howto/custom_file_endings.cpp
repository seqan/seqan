//![includes]
#include <seqan/stream.h>
#include <seqan/seq_io.h>
//![includes]
//![custom_file]
namespace seqan
{
// Your custom file format.
struct MyFastaAdaptor_;
using MyFastaAdaptor = Tag<MyFastaAdaptor_>;

// Specilaize sequence input file with custom tag.
using MySeqFileIn = FormattedFile<Fastq, Input, MyFastaAdaptor>;
//![custom_file]
//![custom_format]
// Your custom format tag.
struct MySeqFormat_;
using MySeqFormat = Tag<MySeqFormat_>;

// The extended TagList containing our custom format.
using MySeqInFormats = TagList<MySeqFormat, SeqInFormats>;

// Overloaded file format metafunction.
template <>
struct FileFormat<FormattedFile<Fastq, Input, MyFastaAdaptor> >
{
    using Type = TagSelector<MySeqInFormats>;
};

// Set magic header.
template <typename T>
struct MagicHeader<MySeqFormat, T> : public MagicHeader<Fasta, T>
{};
//![custom_format]
//![custom_extension]
// Specify the valid ending for your fasta adaptor.
template <typename T>
struct FileExtensions<MySeqFormat, T>
{
    static char const * VALUE[1];
};

template <typename T>
char const * FileExtensions<MySeqFormat, T>::VALUE[1] =
{
    ".fa.dat"  // fasta file with dat ending.
};
//![custom_extension]

//![custom_read_record]
// Overload an inner readRecord function to delegate to the actual fasta parser.
template <typename TIdString, typename TSeqString, typename TSpec>
inline void
readRecord(TIdString & meta, TSeqString & seq, FormattedFile<Fastq, Input, TSpec> & file, MySeqFormat)
{
    readRecord(meta, seq, file.iter, Fasta());  // Just delegate to Fasta parser.
}
} // namespace seqan
//![custom_read_record]
//![main]
int main()
{
    using namespace seqan;
    std::string path = getAbsolutePath("demos/howto/custom_file_ending.fa.dat");

    MySeqFileIn seqFile(path.c_str());

    CharString meta;
    Dna5String seq;

    readRecord(meta, seq, seqFile);

    std::cout << "> " << meta << "\n" << seq << std::endl;
    return 0;
}
//![main]
