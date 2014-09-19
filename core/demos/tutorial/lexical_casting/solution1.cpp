#include <iostream>
#include <seqan/stream.h>

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
            d = seqan::lexicalCast<double>(argv[1]);
            i = seqan::lexicalCast<int>(argv[1]);
            u = seqan::lexicalCast<unsigned>(argv[1]);
        }
        catch (seqan::BadLexicalCast &e)
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
        
        bool bi = seqan::lexicalCast(i, argv[1]);
        bool bu = seqan::lexicalCast(u, argv[1]);
        bool bd = seqan::lexicalCast(d, argv[1]);

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
            seqan::lexicalCastWithException(d, argv[1]);
            seqan::lexicalCastWithException(i, argv[1]);
            seqan::lexicalCastWithException(u, argv[1]);
        }
        catch (seqan::BadLexicalCast &e)
        {
            std::cerr << e.what() << std::endl;
        }

        std::cout << "lexicalCast2<int>(" << argv[1] << ") ==      (" << i << ")\n";
        std::cout << "lexicalCast2<unsigned>(" << argv[1] << ") == (" << u << ")\n";
        std::cout << "lexicalCast2<double>(" << argv[1] << ") ==   (" << d << ")\n";
    }
    
    return 0;
}
