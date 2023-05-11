#include <seqan/sequence.h>

int main(int argc, char const ** argv)
{
    seqan2::String<char> programName = argv[0];
    if (argc > 1)
    {
        seqan2::String<char> firstArg = argv[1];
        if (argc > 2)
            return 1;
    }
    return 0;
}
