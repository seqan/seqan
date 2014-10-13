#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

// Read "<key>\t<value>" map from stdin.  Write out as "<key> -> <value>".

int main()
{
    // Define string and record reader types.  We will read from std::cin which
    // is of type std::istream.  We use a single-pass record reader.
    typedef seqan::DirectionIterator<std::istream, seqan::Input>::Type TReader;

    // Create RecordReader reading from standard input.
    TReader reader = directionIterator(std::cin, seqan::Input());

    // Read the file line by line.
    while (!atEnd(reader))
    {
        // Read first column: The key.
        seqan::CharString key;
        readUntil(key, reader, seqan::IsTab());

        skipOne(reader, seqan::IsTab());    // Skip TAB.

        // Read second column: The value.
        seqan::CharString value;
        readLine(value, reader);    // EOL will not be stored in value.

        // Print ${key} -> ${value}.
        std::cout << key << " -> " << value << std::endl;
    }

    return 0;
}
