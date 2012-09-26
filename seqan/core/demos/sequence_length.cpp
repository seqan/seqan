// File:   sequence_length.cpp
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>

#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/file.h>

using namespace seqan;

int main(int argc, char **argv)
{
  if (argc != 2) {
    std::cerr << "Wrong argument count!" << std::endl
              << "USAGE: sequence_length SEQUENCE.fasta" << std::endl;
    return 1;
  }

  FragmentStore<> fragmentStore;
  if (!loadContigs(fragmentStore, argv[1])) {
    std::cerr << "Could not read contigs." << std::endl;
    return 1;
  }

  size_t total = 0;
  std::cout << "sequence\tlength" << std::endl;
  for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
    std::cout << fragmentStore.contigNameStore[i] << "\t" << length(fragmentStore.contigStore[i].seq) << std::endl;
    total += length(fragmentStore.contigStore[i].seq);
  }
  std::cout << "sum\t" << total << std::endl;

  return 0;
}
