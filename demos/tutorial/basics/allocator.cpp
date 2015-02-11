//![definitions]
#include <iostream>
#include <seqan/basic.h>

using namespace seqan;

int main()
{
    Allocator<MultiPool<> > mpa;
    Allocator<SimpleAlloc<> > sa;
    // store blocksizes in an array
    int bs[3] = {10, 100, 1000};
    int runs = 100000;
    char * buf;
    double startTime;
//![definitions]
//![time-measurements]

    // loop through the different block sizes
    for (int i = 0; i < 3; ++i)
    {
        startTime = sysTime();
        for (int j = 0; j < runs; ++j)
            allocate(mpa, buf, bs[i], TagAllocateTemp());
        clear(mpa);

        std::cout << "Allocating and clearing " << runs << " times blocks of size ";
        std::cout << bs[i] << " with MultiPool Allocator took " << sysTime() - startTime << std::endl;

        startTime = sysTime();
        for (int j = 0; j < runs; ++j)
            allocate(sa, buf, bs[i], TagAllocateTemp());
        clear(sa);

        std::cout << "Allocating and clearing " << runs << " times blocks of size ";
        std::cout << bs[i] << " with Standard Allocator took " << sysTime() - startTime << std::endl;
    }

    return 0;
}
//![time-measurements]
