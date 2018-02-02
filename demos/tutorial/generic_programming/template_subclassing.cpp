// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// This is a simple example of template subclassing.
//
// Referenced by TemplateSubclassing.
// ==========================================================================

#include <iostream>

// Setup class hierarchy by using inheritance for the specialization tags.

struct SpecA {};
struct SpecB {};
struct SpecC :
    SpecB {};
struct SpecD :
    SpecB {};

// Base class -- most generic and thus a fallback.
template <typename TSpec>
struct MyClass
{
    int x;

    MyClass() :
        x(0)
    {}
};

// Specialization for classes B and D.  C automatically inherits everything
// from B and does not overwrite anything.
template <>
struct MyClass<SpecB>
{
    int x;

    MyClass() :
        x(1)
    {}
};

template <>
struct MyClass<SpecD>
{
    int x;

    MyClass() :
        x(2)
    {}
};

// Define some functions that demostrate overloading.

// Most generic case, "base implementation".
template <typename TSpec>
void foo(MyClass<TSpec> const & obj)
{
    std::cout << "foo(MyClass<typename TSpec> const & obj) called!  obj.x == " << obj.x << "\n";
}

// This function overwrites the generic implementation of foo() for the
// specialization B.
template <>
void foo(MyClass<SpecB> const & obj)
{
    std::cout << "foo(MyClass<SpecB> const & obj) called!  obj.x == " << obj.x << "\n";
}

// This function overwrites the generic implementation of foo() for the
// specialization C.
template <>
void foo(MyClass<SpecC> const & obj)
{
    std::cout << "foo(MyClass<SpecC> const & obj) called!  obj.x == " << obj.x << "\n";
}

// The main function calls the functions from above.
int main()
{
    std::cout << "Template Subclassing Demo\n";

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
