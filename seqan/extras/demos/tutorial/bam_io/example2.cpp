#include <seqan/stream.h>

int main()
{
    seqan::Stream<seqan::Bgzf> stream;
    if (!open(stream, "filename.bam", "r"))
        return 1;  // Could not open for reading.
    
    return 0;
}
