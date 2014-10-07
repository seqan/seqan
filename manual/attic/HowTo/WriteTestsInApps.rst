How To: Write Tests in Apps
---------------------------

TOC

This page describes how to write tests for components of apps. These
tests are located within the app folder

Each test program defines a *Test Suite*, a collection of related
*Tests*.

Initial Setup
~~~~~~~~~~~~~

Of course, we first need an app to write the tests for. We create a
simple example sandbox and app, but of course this will also work for
your own app.

::

    #ShellBox
    seqan-trunk # ./util/bin/skel.py repository sandbox/my_sandbox
    ...
    seqan-trunk # ./util/bin/skel.py app my_app sandbox/my_sandbox
    ...

Next, create a header for the code you want to test.

::

    #ShellBox
    seqan-trunk # ./util/bin/skel.py header my_module.h sandbox/my_sandbox/apps/my_app
    ...

Next, add the following function declaration to your header:

::

    #cpp
    // Returns the square of the given parameter.
    int mySquare(int x);

Then, create a file ``sandbox/my_sandbox/apps/my_app/my_module.cpp``
(thus in the same directory as ``module.h``):

::

    #cpp
    #include "my_module.h"

    int mySquare(int x)
    {
        return x * x;  // your code here
    }

To link your application binary against the new module, you have to
replace the ``add_executable()`` call in your
``sandbox/my_sandbox/apps/my_app/CMakeLists.txt``. In the
``CMakeLists.txt``, you should find a line that reads:

::

    add_executable (my_app my_app.cpp)

You have to change this line and add the new module as follows:

::

    add_executable (my_app
                    my_app.cpp
                    my_module.cpp
                    my_module.h)

Next, call CMake (in the build directory) to update your Makefiles or
IDE project files:

::

    #ShellBox
    seqan-build # cmake .
    ...

Creating the Test Program
~~~~~~~~~~~~~~~~~~~~~~~~~

Next, we create the test program and register it in the
``CMakeLists.txt``.

Create a new C++ file at
``sandbox/my_sandbox/apps/my_app/test_my_module.cpp`` with the following
contents:

::

    #cpp

    #undef SEQAN_ENABLE_TESTING
    #define SEQAN_ENABLE_TESTING 1

    #include <seqan/basic.h>

    #include "my_module.h"

    SEQAN_DEFINE_TEST(test_my_module_my_sqr)
    {
        SEQAN_ASSERT_EQ(mySquare(0), 0);
        SEQAN_ASSERT_EQ(mySquare(1), 1);
        SEQAN_ASSERT_EQ(mySquare(2), 4);
        SEQAN_ASSERT_EQ(mySquare(3), 9);
    }

    SEQAN_BEGIN_TESTSUITE(test_suite_name)
    {
        SEQAN_CALL_TEST(test_my_module_my_sqr);
    }
    SEQAN_END_TESTSUITE

Next, create a new executable entry in ``CMakeLists.txt``. We need a new
``add_executable()`` and a new ``target_link_libraries()`` entry for
your new tests.

Add the following below the ``target_link_libraries()`` for ``my_app``:

::

    add_executable (test_my_module
                    test_my_module.cpp
                    my_module.h
                    my_module.cpp)
    target_link_libraries (test_my_module ${SEQAN_LIBRARIES})

Note how similar these lines are to the corresponding lines for
``my_app``.

Next, run CMake in the build directory again to see the
``test_my_module`` program and run it:

::

    #ShellBox
    seqan-build # cmake .
    ...
    seqan-build # make test_my_module && ./bin/test_my_module

Test Suite Skelleton / Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now consider the test program further.

``SEQAN_BEGIN_TESTSUITE(...)`` and ``SEQAN_END_TESTSUITE`` are macros
that expand to book-keeping code for running a test suite.
``SEQAN_DEFINE_TEST(...)`` expands to the definition of a function that
runs a test.

Note that when adding new tests then you have to add them to the
dependencies of the test target in
*sandbox/my\_sandbox/tests/my\_module/CMakeLists.txt*.

Test Macros
~~~~~~~~~~~

Inside your tests, you can use the ``SEQAN_ASSERT*`` and
``SEQAN_ASSERT_*_MSG`` macros to check for assertions. Other useful
macros are seqan:Macro.SEQAN\_PATH\_TO\_ROOT and
seqan:Macro.SEQAN\_TEMP\_FILENAME.

The macros themselves are documented on the `Macro page of the SeqAn
docs <http://www.seqan.de/dddoc/html_devel/INDEXPAGE_Macro.html>`__.

Best Practices
~~~~~~~~~~~~~~

**Rules are there to make you think before you break them.** The
following is not written into stone, but should be good guidelines.
Improvements to the best practices is welcome.

Be Consistent
^^^^^^^^^^^^^

Whatever you do: Be consistent. If the one has read one part of your
code then one should not have to adjust to different variable and
function naming, comment style etc.

Tests Should Compile Without Warnings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    #comment
    TODO(holtgrew): This belongs somewhere else, right?

Make sure that your tests compile without warnings. A common warning is
"comparison of signed and unsigned integer".

In many places, the problematic line looks like this

::

    #cpp
    SEQAN_ASSERT_LT(length(ndl), 30);

The ``length`` function returns an unsigned integer while the string
literal ``30`` represents a (signed) ``int``. You can fix this by
changing the type of the number literal:

::

    #cpp
    SEQAN_ASSERT_LT(length(ndl), 30u);

Break Your Tests Down
^^^^^^^^^^^^^^^^^^^^^

Each test should isolate target an as small as possible and/or feasible
unit of your code. Having short test functions makes them easier to read
and maintain.

Another advantage is that bogus state does not leak into other tests:
Imagine, you have a test that tests a function
``assign_if_positive(a, b)`` that assigns b to a if b is positive.

::

    #cpp
    SEQAN_DEFINE_TEST(test_assign)
    {
        int x = 0;

        assign_if_positive(x, 5);
        SEQAN_ASSERT_EQ(x, 5);

        assign_if_positive(x, -7);
        SEQAN_ASSERT_EQ(x, 5);
    }

Now, what happens if ``assign_if_positive(...)`` has a bug and *never*
assigns a value to its first parameter or always assigns 1? Both of your
assertions will fail. This means you do not really know in which case
the function works well and in which case it does not work well.

Splitting the test make it more robust:

::

    #cpp
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

If you need to initialize the same state for multiple tests, then the
code for this should only exist once. This makes it easier to maintain
since we do not have to change it in multiple places at once. This is
especially useful when following the best practice *Break Your Tests
Down*.

Example:

Instead of

::

    #cpp
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

::

    #cpp

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

Tests can complement examples from the documentation in that they
illustrate each call to your code's API. Thus, make sure that your tests
are well-documented. Not only for users who look up how to use your code
but also for the next maintainer.

There should be a documentation of the test itself and also inline
comments. In your comments, you should focus on the maintainer and not
so much the user. Even if some things are obvious, you might want to
illustrate why you call a function with the given parameters, e.g.
describe the corner cases.

Example:

::

    #cpp

    // Test abs() function with 1, a representant for positive values.
    SEQAN_DEFINE_TEST(test_abs_with_one)
    {
        SEQAN_ASSERT_EQ(abs(1), 1);
    }

    // Test abs() function with 0, the only corner case here.
    SEQAN_DEFINE_TEST(test_abs_with_zero) {
        SEQAN_ASSERT_EQ(abs(0), 0);
    }

    // Test abs() function with -1, a representant for negative values.
    SEQAN_DEFINE_TEST(test_abs_with_minus_one)
    {
        SEQAN_ASSERT_EQ(abs(-1), 1);
    }

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
