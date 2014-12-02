#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

// The following few lines are the actual solution to the assignment.

typedef seqan::OrFunctor<seqan::IsInRange<'0', '9'>,
        seqan::OrFunctor<seqan::IsInRange<'a', 'f'>,
                         seqan::IsInRange<'A', 'F'> > > IsHexDigit;

// This main routine is only some driver code that reads from stdin.

int main()
{
    // We will read from std::cin via an iterator.
    typedef seqan::DirectionIterator<std::istream, seqan::Input>::Type TReader;

    // Create iterator to read from standard input.
    TReader reader = directionIterator(std::cin, seqan::Input());

    seqan::CharString buffer;

    while (!atEnd(reader))
    {
        clear(buffer);
        readUntil(buffer, reader, seqan::NotFunctor<IsHexDigit>());

        // Print hexadecimal number back to the user.
        std::cout << "RECOGNIZED " << buffer << '\n';
        
        // Skip all trailing input.
        skipLine(reader);
    }

    return 0;
}
