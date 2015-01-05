#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

// Read "<key>\t<value>" map from stdin.  Write out as "<key> -> <value>".

int main()
{
    // We will read from std::cin via an iterator.
    typedef DirectionIterator<std::istream, Input>::Type TReader;

    // Create iterator to read from standard input.
    TReader reader = directionIterator(std::cin, Input());

    CharString key, value;

    // Read the file line by line.
    while (!atEnd(reader))
    {
        // Read first column: The key.
        clear(key);
        readUntil(key, reader, IsTab());

        skipOne(reader, IsTab());    // Skip TAB.

        // Read second column: The value.
        clear(value);
        readLine(value, reader);    // EOL will not be stored in value.

        // Print ${key} -> ${value}.
        std::cout << key << " -> " << value << std::endl;
    }

    return 0;
}
