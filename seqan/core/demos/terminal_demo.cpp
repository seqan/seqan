#include <iostream>

#include <seqan/misc/misc_terminal.h>

int main()
{
    unsigned cols = 0, rows = 0;
    seqan::getTerminalSize(cols, rows);
    std::cerr << "cols == " << cols << ", rows == " << rows << std::endl;

    return 0;
}
