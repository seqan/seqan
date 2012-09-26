#include <iostream>
#include <seqan/version.h>

int main()
{
    std::cerr << "SEQAN_VERSION_MAJOR:" << SEQAN_VERSION_MAJOR << "\n"
              << "SEQAN_VERSION_MINOR:" << SEQAN_VERSION_MINOR << "\n"
              << "SEQAN_VERSION_PATCH:" << SEQAN_VERSION_PATCH << "\n"
              << "SEQAN_VERSION_PRE_RELEASE:" << SEQAN_VERSION_PRE_RELEASE << "\n";
    return 0;
}
