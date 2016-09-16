//![assignment1]
#include <seqan/sequence.h>
#include <seqan/index.h>

//![assignment1]

namespace seqan{
//![SAValue]
template<>
struct SAValue<String<Dna> >
{
    typedef unsigned Type;
};
//![SAValue]
//![SAValue2]
template<>
struct SAValue<StringSet<String<Dna> > >
{
    typedef Pair<unsigned, unsigned> Type;
};
//![SAValue2]
//![SAValue3]
typedef Pair<unsigned char, unsigned> Type;
//![SAValue3]
}

//![assignment1]
using namespace seqan;

int main()
{
//![assignment1]
//![save]
    const char * tempFileName = SEQAN_TEMP_FILENAME();
//![save]
    {
//![assignment1]
    String<char> text = "This is the first example";
    Index<String<char>, FMIndex<> > index(text);

//![assignment1]

//![esa]
//![finder]
    String<Dna5> genome = "ACGTACGTACGTN";
    Index<String<Dna5>, IndexEsa<> > esaIndex(genome);
//![esa]
//![finder]
//![fm]
    StringSet<String<AminoAcid> > protein;
    appendValue(protein, "VXLAGZ");
    appendValue(protein, "GKTVXL");
    appendValue(protein, "XLZ");

    Index<StringSet<String<AminoAcid> >, FMIndex<> > fmIndex(protein);
//![fm]
//![finder]
    Finder<Index<String<Dna5>, IndexEsa<> > > esaFinder(esaIndex);
//![finder]
//![finder_multiple]
    find(esaFinder, "ACGT");  // first occurrence of "ACGT"
    find(esaFinder, "ACGT");  // second occurrence of "ACGT"
    find(esaFinder, "ACGT");  // third occurrence of "ACGT"
//![finder_multiple]
//![finder_position]
    find(esaFinder, "ACGT"); // first occurrence of "ACGT"
    position(esaFinder); // -> 0
    find(esaFinder, "ACGT"); // second occurrence of "ACGT"
    position(esaFinder); // -> 4
    find(esaFinder, "ACGT"); // third occurrence of "ACGT"
    position(esaFinder); // -> 8
//![finder_position]

//![save]
    save(index, tempFileName);
//![save]
//![open]
    open(index, tempFileName);
//![open]
//![require]
    indexRequire(index, FibreSA());
//![require]
//![require2]
    indexRequire(esaIndex, EsaSA());
    indexRequire(esaIndex, EsaLcp());
    indexRequire(esaIndex, EsaChildtab());  // for TopDown iterators
    indexRequire(esaIndex, EsaBwt());       // for (Super-)MaxRepeats iterators
//![require2]
//![external]
    Index<String<Dna, External<> >, IndexEsa<> > external_index;
//![external]
    indexRequire(external_index, EsaSA());
    save(external_index, tempFileName);
//![external]
    open(external_index, tempFileName);
//![external]
    StringSet<String<Dna> > set;
    appendValue(set, "ACTG");
//![config]
    typedef FMIndexConfig<void, unsigned> TConfig;
    Index<StringSet<String<Dna> >, FMIndex<void, TConfig> > configIndex(set);
//![config]
    }

    // clean up tmp directory
    ClassTest::_deleteTempFile(ClassTest::_stripFileName(tempFileName));

//![assignment1]
    return 0;
}
//![assignment1]
