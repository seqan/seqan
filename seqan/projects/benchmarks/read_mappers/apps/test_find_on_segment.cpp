#include <seqan/find.h>

using namespace seqan;

// Switch compiler bug on/off.
#define TRIGGER_BUG 1

template <typename THaystack, typename TNeedle>
void doFind(THaystack &myHaystack, TNeedle &myNeedle)
{
  Finder<THaystack> finder(myHaystack);
#if TRIGGER_BUG
  Pattern<TNeedle, DPSearch<Score<int> > > pattern(myNeedle);
#else
  Pattern<TNeedle, MyersUkkonen> pattern(myNeedle);
#endif  // #if TRIGGER_BUG
  setScoreLimit(pattern, -1);
  
  // Find!
  while (find(finder, pattern) and findBegin(finder, pattern)) {
    // nop
    //std::cout << position(finder) << ": " << getScore(pattern) << std::endl;
  }
}

int main(int arg, char **argv)
{
  typedef Dna5String TSequence;
  typedef Segment<TSequence> TSegment;
  typedef ModifiedString<TSequence, ModReverse> TReversedString;
  typedef ModifiedString<TSegment, ModReverse> TReversedModifiedString;
  
  TSequence myHaystack("CGATCGAT");
  TSegment haystackSegment(myHaystack, 1, 6);
  TReversedModifiedString reversedHaystackSegment;
  TSequence myNeedle("CGA");
  TReversedString reversedNeedle(myNeedle);
  
//   doFind(myHaystack, myNeedle);
//   doFind(myHaystack, reversedNeedle);
//   doFind(haystackSegment, myNeedle);
//   doFind(reversedHaystackSegment, myNeedle);
//   doFind(haystackSegment, reversedNeedle);
   doFind(reversedHaystackSegment, reversedNeedle);
}
