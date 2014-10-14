#include <iostream>
#include <fstream>
#include <cstdio>

#include <seqan/stream.h>

// This template function reads the contents from the given Stream in and
// writes it out to std::cout

    std::ifstream in("in.txt", std::ios::binary | std::ios::in);
    std::ofstream out("out.txt", std::ios::binary | std::ios::out);

    // Create iterators to read and write.
    typedef seqan::DirectionIterator<std::ifstream, seqan::Input>::Type TReader;
    typedef seqan::DirectionIterator<std::ofstream, seqan::Output>::Type TWriter;

    TReader reader = directionIterator(in, seqan::Input());
    TWriter writer = directionIterator(out, seqan::Output());

    seqan::CharString buffer;

    while (!atEnd(reader))
    {
        read(buffer, reader, 1000);
        write(writer, buffer);
    }


template <typename TReader>
void doReading(TReader & reader)
{
    seqan::CharString buffer;
    reserve(buffer, 1000);

    while (!atEnd(reader))
    {
        clear(buffer);
        read(buffer, reader, capacity(buffer));
        write(std::cout, buffer);
    }

    return 0;
}

// This template function writes out "Hello World!\n" to the given Stream.

template <typename TWriter>
void doWriting(TWriter & writer)
{
    seqan::CharString buffer = "Hello World!\n";
    write(writer, buffer);
}

// The main function parses the command line, opens the files in the
// appropriate modes with the appropriate stream types and then calls either
// doWriting() or doReading().

int main(int argc, char const ** argv)
{
    if (argc != 4)
    {
        std::cerr << "USAGE: " << argv[0] << " [file|fstream] [r|w] FILENAME\n";
        return 1;
    }

    // Check first argument.
    if (seqan::CharString(argv[1]) != "file" && seqan::CharString(argv[1]) != "fstream")
    {
        std::cerr << "ERROR: " << argv[1] << " is not a valid stream type name.\n";
        return 1;
    }
    bool useFile = (seqan::CharString(argv[1]) == "file");

    // Check second argument.
    if (seqan::CharString(argv[2]) != "r" && seqan::CharString(argv[2]) != "w")
    {
        std::cerr << "ERROR: " << argv[2] << " is not a valid operation name.\n";
        return 1;
    }
    bool doRead = (seqan::CharString(argv[2]) == "r");

    // Branches for stream and operation type.
    int res = 0;
    if (useFile)  // FILE *
    {
        FILE * fp;
        
        if (doRead)  // reading
            fp = fopen(argv[3], "rb");
        else  // writing
            fp = fopen(argv[3], "wb");

        if (fp == 0)
        {
            std::cerr << "ERROR: Could not open " << argv[3] << "\n";
            return 1;
        }

        if (doRead)  // reading
            res = doReading(fp);
        else  // writing
            res = doWriting(fp);

        fclose(fp);
    }
    else  // std::fstream
    {
        std::fstream stream;
        
        if (doRead)  // reading
            stream.open(argv[3], std::ios::binary | std::ios::in);
        else  // writing
            stream.open(argv[3], std::ios::binary | std::ios::out);

        if (!stream.good())
        {
            std::cerr << "ERROR: Could not open " << argv[3] << "\n";
            return 1;
        }

        if (doRead)  // reading
            res = doReading(stream);
        else  // writing
            res = doWriting(stream);
    }

    if (res != 0)
        std::cerr << "ERROR: There was an error accessing the file!\n";
    return res;
}
