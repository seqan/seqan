#include <seqan/stream.h>

int main()
{
    seqan::Stream<seqan::Bgzf> stream;
    if (!open(stream, "filename.bam", "w"))
        return 1;  // Could not open for writing.
    
    return 0;
}
