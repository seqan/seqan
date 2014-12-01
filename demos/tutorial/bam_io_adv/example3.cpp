#include <seqan/bam_io.h>

int main()
{
    seqan::BamFileOut file;
    if (!open(file, "filename.bam"))
        return 1;  // Could not open for writing.
    
    return 0;
}
