#include <iostream>

#include <seqan/misc/misc_terminal.h>

int main()
{
    // Get terminal size and print it to stdout.
    unsigned cols = 0, rows = 0;
    seqan::getTerminalSize(cols, rows);
    std::cout << "cols == " << cols << ", rows == " << rows << "\n";

    return 0;
}
