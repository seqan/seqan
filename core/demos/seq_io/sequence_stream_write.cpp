#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char ** argv)
{
    CharString path = SEQAN_TEMP_FILENAME();
    append(path, ".fa");

    SequenceFile<Output> file(toCString(path));

    CharString meta = "chr1";
    CharString seq = "ACGT";

    write(file, meta, seq);

    close(file);

    return 0;
}
