//![header]
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>

#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;
//![header]

//![main]
int main(int argc, char const ** argv)
{
    typedef VirtualStream<char, Input> TVStream;

    if (argc != 2)
    {
        std::vector<std::string> file_ext_list = TVStream::getFileExtensions();
        std::ostream_iterator<std::string> out_it (std::cerr,"|");

        std::cerr << "USAGE: " << argv[0] << " input[";
        std::copy(file_ext_list.begin(), file_ext_list.end() - 1, out_it);
        std::cerr << "]" << std::endl;
        return 1;
    }

    TVStream vin;

    if (!open(vin, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << std::endl;
        return 1;
    }

    // Create iterators to read and write.
    typedef DirectionIterator<TVStream, Input>::Type TReader;
    TReader reader = directionIterator(vin, Input());

    CharString buffer;
    reserve(buffer, 1000);

    while (!atEnd(reader))
    {
        clear(buffer);
        read(buffer, reader, capacity(buffer));
        write(std::cout, buffer);
    }

    return 0;
}
//![main]
