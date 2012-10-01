/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2007
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  This is a simple example of template subclassing.
  ===========================================================================*/

#include <iostream>

// Setup class hierarchy by using inheritance for the specialization tags.
struct SpecA {};
struct SpecB {};
struct SpecC : SpecB {};
struct SpecD : SpecB {};


// Base class -- most generic and thus a fallback.
template <typename TSpec>
struct MyClass {
    int x;

    MyClass() : x(0) {}
};


// Specialization for classes B and D.  C automatically inherits
// everything from B and does not overwrite anything.
template <>
struct MyClass<SpecB> {
    int x;

    MyClass() : x(1) {}
};


template <>
struct MyClass<SpecD> {
    int x;

    MyClass() : x(2) {}
};


// Define some functions that demostrate overloading.

// Most generic case, "base implementation".
template <typename TSpec>
void foo(MyClass<TSpec> const & obj) {
    std::cout << "foo(MyClass<typename TSpec> const & obj) called!  obj.x == " << obj.x << std::endl;
}


// This function overwrites the generic implementation of foo() for
// the specialization B.
template <>
void foo(MyClass<SpecB> const & obj) {
    std::cout << "foo(MyClass<SpecB> const & obj) called!  obj.x == " << obj.x << std::endl;
}


// This function overwrites the generic implementation of foo() for
// the specialization C.
template <>
void foo(MyClass<SpecC> const & obj) {
    std::cout << "foo(MyClass<SpecC> const & obj) called!  obj.x == " << obj.x << std::endl;
}


int main(int, char**) {
    std::cout << "Template Subclassing Demo" << std::endl;

    MyClass<SpecA> a;
    MyClass<SpecB> b;
    MyClass<SpecC> c;
    MyClass<SpecD> d;

    foo(a);  // => foo(MyClass<typename TSpec> const & obj) called!  obj.x == 0
    foo(b);  // => foo(MyClass<SpecB> const & obj) called!  obj.x == 1
    foo(c);  // => foo(MyClass<SpecC> const & obj) called!  obj.x == 0
    foo(d);  // => foo(MyClass<typename TSpec> const & obj) called!  obj.x == 2
    
    return 0;
}
