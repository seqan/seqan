#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

// The following few lines are the actual solution to the assignment.

struct Punct_;
typedef seqan::Tag<Punct_> Punct;

inline int
_charCompare(int const c, Punct const & /* tag*/)
{
    return ispunct(c);
}

template <typename TStream, typename TPass, typename TBuffer>
inline int
readPunctuation(TBuffer & buffer, seqan::RecordReader<TStream, TPass> & reader)
{
    return seqan::_readHelper(buffer, reader, Punct(), false);
}

// This main routine is only some driver code that reads from stdin.

int main()
{
    seqan::RecordReader<std::istream, seqan::SinglePass<> > reader(std::cin);

    while (!atEnd(reader))
    {
        seqan::CharString buffer;
        int res = readPunctuation(buffer, reader);
        if (res != 0 && res != seqan::EOF_BEFORE_SUCCESS)
        {
            std::cerr << "ERROR: Could not read from standard input.\n";
            return 1;
        }

        // Print hexadecimal number back to the user.
        std::cout << "RECOGNIZED " << buffer << '\n';
        
        // Skip all trailing input.
        skipLine(reader);
    }

    return 0;
}
