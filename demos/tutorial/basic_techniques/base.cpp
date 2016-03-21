#include <seqan/basic.h>
//![class]
template <typename T>
class vector{};
//![class]

using namespace seqan;

int main()
{
//![vector]
    vector<int> my_vector;
//![vector]
    ignoreUnusedVariableWarning(my_vector);
    return 0;
}
