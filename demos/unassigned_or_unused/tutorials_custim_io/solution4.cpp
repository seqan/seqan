#include <iostream>
#include <seqan/stream.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: tutorial_solution1 VALUE\n";
        return 1;
    }

    // Lexical casting with the 1-argument lexicalCast().
    {
        int      i = 0;
        unsigned u = 0;
        double   d = 0;

        try
        {
            d = lexicalCast<double>(argv[1]);
            i = lexicalCast<int>(argv[1]);
            u = lexicalCast<unsigned>(argv[1]);
        }
        catch (BadLexicalCast & e)
        {
            std::cerr << e.what() << std::endl;
        }
        std::cout << "lexicalCast<int>(" << argv[1] << ") ==      " << i << '\n';
        std::cout << "lexicalCast<unsinged>(" << argv[1] << ") == " << u << '\n';
        std::cout << "lexicalCast<double>(" << argv[1] << ") ==   " << d << '\n';
    }

    // Lexical casting with the 2-argument lexicalCast().
    {
        int      i = 0;
        unsigned u = 0;
        double   d = 0;

        bool bi = lexicalCast(i, argv[1]);
        bool bu = lexicalCast(u, argv[1]);
        bool bd = lexicalCast(d, argv[1]);

        std::cout << "lexicalCast2<int>(" << argv[1] << ") ==      (" << bi << ", " << i << ")\n";
        std::cout << "lexicalCast2<unsigned>(" << argv[1] << ") == (" << bu << ", " << u << ")\n";
        std::cout << "lexicalCast2<double>(" << argv[1] << ") ==   (" << bd << ", " << d << ")\n";
    }

    // Lexical casting with the 2-argument lexicalCast() that throws exceptions.
    {
        int      i = 0;
        unsigned u = 0;
        double   d = 0;

        try
        {
            lexicalCastWithException(d, argv[1]);
            lexicalCastWithException(i, argv[1]);
            lexicalCastWithException(u, argv[1]);
        }
        catch (BadLexicalCast & e)
        {
            std::cerr << e.what() << std::endl;
        }

        std::cout << "lexicalCast2<int>(" << argv[1] << ") ==      (" << i << ")\n";
        std::cout << "lexicalCast2<unsigned>(" << argv[1] << ") == (" << u << ")\n";
        std::cout << "lexicalCast2<double>(" << argv[1] << ") ==   (" << d << ")\n";
    }

    return 0;
}
