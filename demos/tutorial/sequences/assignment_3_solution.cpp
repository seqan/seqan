   #include <seqan/stream.h>
   #include <seqan/sequence.h>
   #include <seqan/file.h>

   using namespace seqan;

   int main()
   {
       String<Dna5> nucleotides = "AGTCGTGNNANCT";
       String<Dna5> lesser;
       String<Dna5> greater;

       for (unsigned i = 0; i < length(nucleotides); ++i){
	   if (nucleotides[i] < 'G')
	       appendValue(lesser, nucleotides[i]);
	   else if (nucleotides[i] > 'G')
	       appendValue(greater, nucleotides[i]);
       }
       std::cout << "Lesser nucleotides: " << lesser << std::endl;
       std::cout << "Greater nucleotides: " << greater << std::endl;
   }

