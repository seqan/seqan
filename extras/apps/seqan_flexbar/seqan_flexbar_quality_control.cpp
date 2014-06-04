#include "seqan_flexbar.h"

int main(int argc, char const ** argv)
{
    // Run quality control program.
    flexiProgram = QUALITY_CONTROL;
    return flexbarMain(argc, argv);
}
