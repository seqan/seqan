#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
   //![assign]
   String<Dna> dna_source = "acgtgcat";
   String<char> char_target;
   assign(char_target, dna_source);
   std::cout << "Copy Conversion: target sequence after assignment: " << char_target << std::endl;

   clear(char_target);
   std::cout << dna_source << std::endl;
   move(char_target, dna_source);
   std::cout << "Copy Conversion: target sequence after assignment: " << char_target << std::endl;
   std::cout << dna_source << std::endl;
   //![assign]

   //![move]
   String<char> char_source = "acgtgcat";
   String<Dna> dna_target;

   // The in-place move conversion.
   move(dna_target, char_source);
   std::cout << "Move Conversion: target sequence after moving: " << dna_target << std::endl;
   //![move]
}
