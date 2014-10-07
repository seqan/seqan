#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    CharString path = SEQAN_TEMP_FILENAME();
    append(path, ".fa");

    SeqFileOut file(toCString(path));
    writeRecord(file, "chr1", "ACGT");
    close(file);

    return 0;
}
