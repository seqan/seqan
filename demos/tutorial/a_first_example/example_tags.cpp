template <typename T>
struct Tag{};

struct TagA_;
typedef Tag<TagA_> TagA;

struct TagB_;
typedef Tag<TagB_> TagB;

void myFunction(TagA const &){}  // (1)
void myFunction(TagB const &){}  // (2)

int main()
{
    myFunction(TagA());  // (3)
    myFunction(TagB());  // (4)
    return 0;
}
