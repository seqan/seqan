//![main]
#include <seqan/stream.h>
#include <seqan/sequence_journaled.h>

using namespace seqan;

int main()
{
//![main]
//![typedef]
    typedef String<char, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournaledString;
    typedef Host<TJournaledString>::Type THost;
//![typedef]

//![init]
    THost hostStr = "thisisahostsequence";
    TJournaledString journalStr;
    setHost(journalStr, hostStr);

    std::cout << "After creating the Journaled String:" << std::endl;
    std::cout << "Host: " << host(journalStr) << std::endl;
    std::cout << "Journal: " << journalStr << std::endl;
    std::cout << "Nodes: " << journalStr._journalEntries << std::endl;
    std::cout << std::endl;
//![init]

//![modification]
    insert(journalStr, 7, "modified");
    erase(journalStr, 19, 27);

    std::cout << "After modifying the Journaled String:" << std::endl;
    std::cout << "Host: " << host(journalStr) << std::endl;
    std::cout << "Journal: " << journalStr << std::endl;
    std::cout << "Nodes: " << journalStr._journalEntries << std::endl;
    std::cout << std::endl;
//![modification]

//![flatten]
    flatten(journalStr);
    std::cout << "After flatten the Journaled String:" << std::endl;
    std::cout << "Host: " << host(journalStr) << std::endl;
    std::cout << "Journal: " << journalStr << std::endl;
    std::cout << "Nodes: " << journalStr._journalEntries << std::endl;

    return 0;
}
//![flatten]
