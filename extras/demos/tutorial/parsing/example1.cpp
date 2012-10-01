#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

// Read "<key>\t<value>" map from stdin.  Write out as "<key> -> <value>".

int main()
{
    // Define string and record reader types.  We will read from std::cin which
    // is of type std::istream.  We use a single-pass record reader.
    typedef seqan::RecordReader<std::istream, seqan::SinglePass<> > TRecordReader;

    int res = 0;  // Used to store I/O results.

    // Create RecordReader reading from standard input.
    TRecordReader reader(std::cin);

    // Read the file line by line.
    while (!atEnd(reader))
    {
        // Read first column: The key.
        seqan::CharString key;
        res = readUntilChar(key, reader, '\t');
        if (res != 0)
            return 1;

        goNext(reader);  // Skip TAB.

        // Read second column: The value.
        seqan::CharString value;
        res = readLine(value, reader);  // EOL will not be stored in value.
        if (res != 0)
            return 1;

        // Print ${key} -> ${value}.
        std::cout << key << " -> " << value << std::endl;
    }

    return 0;
}
