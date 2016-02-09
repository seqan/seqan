struct SpecA;
struct SpecB;
struct SpecC;

template <typename TAlphabet, typename TSpec>
class String{};

template <typename TAlphabet, typename TSpec>
void myFunction(String<TAlphabet, TSpec> const &){}  // Variant (A)

template <typename TAlphabet>
void myFunction(String<TAlphabet, SpecB> const &){}  // Variant (B)

// ...

int main()
{
    String<char, SpecA> a;
    String<char, SpecB> b;
    String<char, SpecC> c;

    myFunction(a);            // calls (A)
    myFunction(b);            // calls (B)
    myFunction(c);            // calls (A)
}
