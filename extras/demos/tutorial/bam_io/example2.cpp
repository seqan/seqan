#include <seqan/bam_io.h>

int main()
{
    seqan::BamFileIn file;
    if (!open(file, "filename.bam"))
        return 1;  // Could not open for reading.
    
    return 0;
}
