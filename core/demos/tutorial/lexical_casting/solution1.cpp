#include <iostream>
#include <seqan/stream.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: tutorial_solution1 VALUE\n";
        return 1;
    }

    // Lexical casting with lexicalCast().
    {
        int      i = 0;
        unsigned u = 0;
        double   d = 0;
        
        i = seqan::lexicalCast<int>(argv[1]);
        u = seqan::lexicalCast<int>(argv[1]);
        d = seqan::lexicalCast<int>(argv[1]);
        std::cout << "lexicalCast<int>(" << argv[1] << ") ==      " << i << '\n';
        std::cout << "lexicalCast<unsinged>(" << argv[1] << ") == " << u << '\n';
        std::cout << "lexicalCast<double>(" << argv[1] << ") ==   " << d << '\n';
    }
    
    // Lexical casting with lexicalCast2().
    {
        int      i = 0;
        unsigned u = 0;
        double   d = 0;
        
        bool bi = seqan::lexicalCast2(i, argv[1]);
        bool bu = seqan::lexicalCast2(u, argv[1]);
        bool bd = seqan::lexicalCast2(d, argv[1]);

        std::cout << "lexicalCast2<int>(" << argv[1] << ") ==      (" << bi << ", " << i << ")\n";
        std::cout << "lexicalCast2<unsigned>(" << argv[1] << ") == (" << bu << ", " << u << ")\n";
        std::cout << "lexicalCast2<double>(" << argv[1] << ") ==   (" << bd << ", " << d << ")\n";
    }
    
    return 0;
}
