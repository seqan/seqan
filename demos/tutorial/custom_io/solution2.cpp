#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

// The following few lines are the actual solution to the assignment.

typedef OrFunctor<IsInRange<'0', '9'>,
                  OrFunctor<IsInRange<'a', 'f'>,
                            IsInRange<'A', 'F'> > > IsHexDigit;

// This main routine is only some driver code that reads from stdin.

int main()
{
    // We will read from std::cin via an iterator.
    typedef DirectionIterator<std::istream, Input>::Type TReader;

    // Create iterator to read from standard input.
    TReader reader = directionIterator(std::cin, Input());

    CharString buffer;

    while (!atEnd(reader))
    {
        clear(buffer);
        readUntil(buffer, reader, NotFunctor<IsHexDigit>());

        // Print hexadecimal number back to the user.
        std::cout << "RECOGNIZED " << buffer << '\n';

        // Skip all trailing input.
        skipLine(reader);
    }

    return 0;
}
