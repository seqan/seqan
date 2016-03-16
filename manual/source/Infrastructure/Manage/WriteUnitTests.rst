.. sidebar:: ToC

    .. contents::

.. _infra-manage-write-unit-tests:

Writing Unit Tests
==================

This page describes how to write tests for the SeqAn library.
Each test program defines a *Test Suite*, a collection of related *Tests*.

Test Suite Skeleton / Example
------------------------------

A skeleton and example for a test suite program looks as follows:

.. code-block:: cpp

   #include <seqan/basic.h>

   SEQAN_DEFINE_TEST(test_suite_name_test_name)
   {
       int ii = 1;
       for (int jj = 0; jj < 10; ++jj)
       {
           ii *= 2;
       }
       SEQAN_ASSERT_EQ(ii, 1024);
   }

   SEQAN_BEGIN_TESTSUITE(test_suite_name)
   {
       SEQAN_CALL_TEST(test_suite_name_test_name);
   }
   SEQAN_END_TESTSUITE

``SEQAN_BEGIN_TESTSUITE(...)`` and ``SEQAN_END_TESTSUITE`` are macros that expand to book-keeping code for running a test suite.
``SEQAN_DEFINE_TEST(...)`` expands to the definition of a function that runs a test.

Getting Started With Our Test Template
--------------------------------------

To make creating tests easier the code generator ``util/bin/skel.py`` has a command to generate test skeletons for you.
As parameters, you give it the name of the module you want to test and the path to the repository.
For example, use ``skel.py tests my_module .`` to create tests for the module *my_module* in the directory ``tests``:

.. code-block:: console

   seqan $ ./util/bin/skel.py test my_module .
   ...
   tests/my_module/
   ├── CMakeLists.txt
   ├── test_my_module.cpp
   └── test_my_module.h

Afterwards, you can compile and run the tests:

.. code-block:: console

   $ mkdir -p build/Debug
   $ cd build/Debug
   $ cmake ../..
   $ make test_my_module
   $ ./tests/my_module/test_my_module
   ...

.. note::

   When adding new tests you have to add them to the dependencies of the test target in *tests/my_module/CMakeLists.txt*.

Test Macros
-----------

Inside your tests, you can use the ``SEQAN_ASSERT*`` and ``SEQAN_ASSERT_*_MSG`` macros to check for assertions.
Other useful macros are :dox:`SEQAN_PATH_TO_ROOT` and :dox:`SEQAN_TEMP_FILENAME`.

The macros themselves are documented in the dox: ``SeqAn API documentation AssertMacros``.

Assertion Caveats
-----------------

When using one of the LT/GT/LEQ/GEQ/EQ/NEQ macros, the values have to provide a stream operator (``operator<<``) to write them to an output stream.
If this is not implemented, then the assertion will not compile and something like the following will be printed by the compiler (in this case the GCC).

.. code-block:: console

   In file included from seqan/basic.h:55:0,
                    from tests/sequence/test_sequence.cpp:4:
   seqan/basic/basic_testing.h: In function 'bool ClassTest::testEqual(const char*, int, const T1&, const char*, const T2&, const char*, const char*, ...) [with T1 = Iter<String<char, Block<3u> >, PositionIterator>, T2 = Iter<String<char, Block<3u> >, PositionIterator>]':
   seqan/basic/basic_testing.h:435:81:   instantiated from 'bool ClassTest::testEqual(const char*, int, const T1&, const char*, const T2&, const char*) [with T1 = Iter<String<char, Block<3u> >, PositionIterator>, T2 = Iter<String<char, Block<3u> >, PositionIterator>]'
   tests/sequence/test_string.h:386:2:   instantiated from 'void TestStringBasics() [with TMe = String<char, Block<3u> >]'
   tests/sequence/test_string.h:475:45:   instantiated from here
   seqan/basic/basic_testing.h:385:13: error: no match for 'operator<<' in 'std::operator<< [with _Traits = std::char_traits<char>](((std::ostream&)((std::ostream*)std::operator<< [with _Traits = std::char_traits<char>](((std::ostream&)((std::ostream*)std::operator<< [with _Traits = std::char_traits<char>](((std::ostream&)((std::ostream*)std::operator<< [with _Traits = std::char_traits<char>](((std::ostream&)((std::ostream*)std::operator<< [with _Traits = std::char_traits<char>](((std::ostream&)((std::ostream*)((std::ostream*)std::operator<< [with _Traits = std::char_traits<char>](((std::ostream&)((std::ostream*)std::operator<< [with _Traits = std::char_traits<char>](((std::ostream&)(& std::cerr)), file))), ((const char*)":")))->std::basic_ostream<_CharT, _Traits>::operator<< [with _CharT = char, _Traits = std::char_traits<char>](line))), ((const char*)" Assertion failed : ")))), expression1))), ((const char*)" == ")))), expression2))), ((const char*)" was: ")) << value1'

The workaround is to use

.. code-block:: cpp

   SEQAN_ASSERT(end(str3) == begin(str3) + 7);

instead of

.. code-block:: cpp

    SEQAN_ASSERT_EQ(end(str3), begin(str3) + 7);

Best Practices
--------------

**Rules are there to make you think before you break them.**
The following is not written into stone, but should be good guidelines.
Improvements to the best practices is welcome.

Be Consistent
^^^^^^^^^^^^^

Whatever you do: Be consistent.
If the one has read one part of your code then one should not have to adjust to different variable and function naming, comment style etc.

Tests Should Compile Without Warnings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure that your tests compile without warnings.
A common warning is "comparison of signed and unsigned integer".

In many places, the problematic line looks like this

.. code-block:: cpp

   SEQAN_ASSERT_LT(length(ndl), 30);

The ``length`` function returns an unsigned integer while the string literal ``30`` represents a (signed) ``int``.
You can fix this by changing the type of the number literal:

.. code-block:: cpp

    SEQAN_ASSERT_LT(length(ndl), 30u);

Break Your Tests Down
^^^^^^^^^^^^^^^^^^^^^

Each test should verify a part of the library as small as possible while still being meaningful.
Having short test functions makes them easier to read and maintain.

Another advantage is that bogus state does not leak into other tests: imagine, you have a test that tests a function ``assign_if_positive(a, b)`` that assigns b to a if b is positive.

.. code-block:: cpp

   SEQAN_DEFINE_TEST(test_assign)
   {
       int x = 0;

       assign_if_positive(x, 5);
       SEQAN_ASSERT_EQ(x, 5);

       assign_if_positive(x, -7);
       SEQAN_ASSERT_EQ(x, 5);
   }

Now, what happens if ``assign_if_positive(...)`` has a bug and *never* assigns a value to its first parameter or always assigns 1?
Both of your assertions will fail.
This means you do not really know in which case the function works well and in which case it does not work well.

Splitting the test makes it more robust:

.. code-block:: cpp

   SEQAN_DEFINE_TEST(test_assign_positive)
   {
       int x = 0;
       assign_if_positive(x, 5);
       SEQAN_ASSERT_EQ(x, 5);
   }

   SEQAN_DEFINE_TEST(test_assign_negative)
   {
       int x = 0;
       assign_if_positive(x, -7);
       SEQAN_ASSERT_EQ(x, 0);
   }

Use Helper Functions For Setup/TearDown
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you need to initialize the same state for multiple tests, then the code for this should only exist once.
This makes it easier to maintain since we do not have to change it in multiple places at once.
This is especially useful when following the best practice `Break Your Tests Down`_.

Example:

Instead of

.. code-block:: cpp

   SEQAN_DEFINE_TEST(test_grep)
   {
       char *contents = loadFile("corpus.txt");

       int pos = doGrep(contents, "nonexisting pattern");
       SEQAN_ASSERT_EQ(pos, -1);

       pos = doGrep(contents, "existing pattern");
       SEQAN_ASSERT_EQ(pos, 3);

       delete contents;
   }

do

.. code-block:: cpp

   // Set-up for test_grep_{success, failure}.
   void testGrepSetUp(const char *filename, char *outContents)
   {
       outContents = loadFile(filename);
   }

   // Tear-down for test_grep_{success, failure}.
   void testGraphTearDown(char *contents)
   {
       delete contents;
   }

   // Test greping for existing patterns.
   SEQAN_DEFINE_TEST(test_grep_success)
   {
       // corpus.txt contains the string "1234existing pattern567".
       char *contents;
       testGrepSetUp("corpus.txt", contents);

       int pos = doGrep(contents, "existing pattern");
       SEQAN_ASSERT_EQ(pos, 3);

       testGrepTearDown(contents);
   }

   // Test greping for non-existing patterns.
   SEQAN_DEFINE_TEST(test_grep_failure)
   {
       // corpus.txt contains the string "1234existing pattern567".
       char *contents;
       testGrepSetUp("corpus.txt", contents);

       int pos = doGrep(contents, "nonexisting pattern");
       SEQAN_ASSERT_EQ(pos, -1);

       testGrepTearDown(contents);
   }

Comment Your Tests
^^^^^^^^^^^^^^^^^^

Tests can complement examples from the documentation in that they illustrate each call to your code's API.
Thus, make sure that your tests are well-documented.
Not only for users who look up how to use your code but also for the next maintainer.

There should be a documentation of the test itself and also inline comments.
In your comments, you should focus on the maintainer and not so much on the user.
Even if some things are obvious, you might want to illustrate why you call a function with the given parameters, e.g. describe the corner cases.

Example:

.. code-block:: cpp

   // Test abs() function with 1, a representative for positive values.
   SEQAN_DEFINE_TEST(test_abs_with_one)
   {
       SEQAN_ASSERT_EQ(abs(1), 1);
   }

   // Test abs() function with 0, the only corner case here.
   SEQAN_DEFINE_TEST(test_abs_with_zero)
   {
       SEQAN_ASSERT_EQ(abs(0), 0);
   }

   // Test abs() function with -1, a representative for negative values.
   SEQAN_DEFINE_TEST(test_abs_with_minus_one)
   {
       SEQAN_ASSERT_EQ(abs(-1), 1);
   }
