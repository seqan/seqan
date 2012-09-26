#include <iostream>

#include <seqan/stream.h>

int main(int, char const **)
{
    using namespace seqan;

    int resultI = 0;
    double resultD = 0;
    bool b = false;

    resultI = lexicalCast<int>("123");
    std::cerr << "lexicalCast<int>(\"123\")   --> " << resultI << std::endl;

    resultI = lexicalCast<int>("123XX");
    std::cerr << "lexicalCast<int>(\"123XX\") --> " << resultI << std::endl;

    b = lexicalCast2<int>(resultI, "-123");
    std::cerr << "lexicalCast2<int>(\"-123\") --> (" << b << ", " << resultI << ")" << std::endl;
    
    b = lexicalCast2<double>(resultD, "-123");
    std::cerr << "lexicalCast2<double>(\"-123\") --> (" << b << ", " << resultD << ")" << std::endl;
    
    return 0;
}
