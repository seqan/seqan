#include "seqan_flexbar.h"

int main(int argc, char const ** argv)
{
    // Run quality control program.
    flexiProgram = ADAPTER_REMOVAL;
    return flexbarMain(argc, argv);
}
