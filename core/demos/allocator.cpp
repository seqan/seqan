///A tutorial about the use of allocators.

#include <seqan/basic.h>
using namespace seqan;

///We define an arbitrary class.
struct MyClass
{
};

int main()
{
///We create 100 instances of $MyClass$ on the heap 
///using a default temporary allocator object $Default$.
	MyClass* my_class_arr;
	allocate(Default(), my_class_arr, 100);
	arrayConstruct(my_class_arr, my_class_arr + 100);
///Before the storage is deallocated, the $MyClass$ objects have to be destroyed.
	arrayDestruct(my_class_arr, my_class_arr + 100);
	deallocate(Default(), my_class_arr, 100);
///We can use any kind of object as an allocator.
///However, dedicated allocators offer more advanced functionality, e.g. @Function.clear@.
	Allocator<SimpleAlloc< > > alloc1;
	char * char_array;
	allocate(alloc1, char_array, 300);
///@Function.clear@ can be used to deallocate all storage at once.	
	clear(alloc1);
	return 0;
}
