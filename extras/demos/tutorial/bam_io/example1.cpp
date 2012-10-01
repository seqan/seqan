#include <seqan/bam_io.h>
#include <seqan/sequence.h>

int main()
{
    // Create some shortcuts to the types that we will use.
    typedef seqan::StringSet<seqan::CharString> TNameStore;
    typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
    typedef seqan::BamIOContext<TNameStore>     TBamIOContext;

    // Setup the variables.
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   bamIOContext(nameStore, nameStoreCache);

    return 0;
}
