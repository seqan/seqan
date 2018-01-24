.. sidebar:: ToC

    .. contents::

.. _tutorial-getting-started-first-steps-in-seqan:

A First Example
===============

Learning Objective
  You will learn the most basic concepts of SeqAn.
  After this tutorial you will be ready to deal with the more specific tutorials, e.g. Sequences.

Difficulty
  Very basic

Duration
  1.5h

Prerequisites
  Basic C or C++ knowledge

Welcome to the SeqAn "Hello World".
This is the first practical tutorial you should look at when starting to use our software library.

We assume that you have some programming experience (preferably in C++ or C) and concentrate on SeqAn specific aspects.
We will start out pretty slowly and hopefully the tutorial will make sense to you even if you are new to C++.
However, to really leverage the power of SeqAn you will have to learn C++.
There are many tutorials on C++, for example `the tutorial at cplusplus.com <http://www.cplusplus.com/doc/tutorial/>`_.

This tutorial will walk you through a simple example program that highlights the things that are most prominently different from the libraries that many SeqAn newcomers are used to:

* extensive usage of C++ templates,
* generic programming using templates,
* using references instead of pointers in most places,
* and more.

Running Example
---------------

Let's start with a simple example programm. The program will do a pattern search of a short query sequence (pattern) in a long subject sequence (text).
We define the score for each position of the database sequence as the sum of matching characters between the pattern and the text.

The following figure shows an expected result:

::

    score:    101 ...        ... 801 ...
    text:     This is an awesome tutorial to get to know SeqAn!
    pattern:  tutorial           tutorial
               tutorial           tutorial
                ...                ...


The first position has a score of 1, because the ``i`` in the pattern matches the ``i`` in ``is``.
This is only a toy example for explanatory reasons and we ignore any more advanced implementations.

In SeqAn the program could look like this (we will explain every line of code shortly):

.. includefrags:: demos/tutorial/a_first_example/basic_code.cpp
   :fragment: all

Whenever we use SeqAn classes or functions we have to explicitly write the namespace qualifier ``seqan::`` in front of the class name or function.
This can be circumvented if we include the line ``using namespace seqan;`` at the top of the working example.
However, during this tutorial we will not do this, such that SeqAn classes and functions can be recognized more easily.

.. attention::

   Argument-Dependent Name Lookup (Koenig Lookup)

   Using the namespace prefix ``seqan::`` is not really necessary in all places.
   In many cases, the Koenig lookup rule in C++ for functions makes this unnecessary.
   Consider the following, compiling, example.

   .. includefrags:: demos/tutorial/a_first_example/base.cpp
       :fragment: lookup_rule

   Here, the function ``length`` does not have a namespace prefix.
   The code compiles nevertheless.
   The compiler automatically looks for a function ``length`` in the namespace of its arguments.

Note that we follow the rules for variable, function, and class names as outlined in the :ref:`SeqAn style guide <infra-contribute-style-cpp>`.
For example:
1. variables and functions use lower case,
2. struct, enum and classes use CamelCase,
3. metafunctions start with a capital letter, and
4. metafunction values are UPPERCASE.

Assignment 1
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Create a demo program and replace its content with the code above.

   Hint
     Depending on your operating system you have different alternatives to create a demo application.
     An in depth description can be found in GettingStarted.

   Solution
     Click ''more...''

     .. container:: foldable

        .. includefrags:: demos/tutorial/a_first_example/solution_1.cpp

        Your output should look like this:

        .. includefrags:: demos/tutorial/a_first_example/solution_1.cpp.stdout

SeqAn and Templates
-------------------

Let us now have a detailed look at the program.

We first include the IOStreams library that we need to print to the screen and the SeqAn's ``<seqan/file.h>`` as well as ``<seqan/sequence.h>`` module from the SeqAn library that provides SeqAn :dox:`String`.

.. includefrags:: demos/tutorial/a_first_example/basic_code_detailed.cpp
   :fragment: includes

The :dox:`String String class` is one of the most fundamental classes in SeqAn, which comes as no surprise since SeqAn is used to analyse sequences (there is an extra tutorial for SeqAn :ref:`sequences <tutorial-datastructures-sequences>` and :ref:`alphabets <tutorial-datastructures-sequences-alphabets>`).

In contrast to the popular string classes of Java or C++, SeqAn provides different string implementations and different alphabets for its strings.
There is one string implementation that stores characters in memory, just like normal C++ strings.
Another string implementation stores the characters on disk and only keeps a part of the sequence in memory.
For alphabets, you can use strings of nucleotides, such as genomes, or you can use strings of amino acids, for example.

SeqAn uses **template functions** and **template classes** to implement the different types of strings using the **generic programming** paradigm.
Template functions/classes are normal functions/classes with the additional feature that one passes the type of a variable as well as its value (see also: `templates in cpp <http://www.cplusplus.com/doc/tutorial/templates/>`_).
This means that SeqAn algorithms and data structures are implemented in such a way that they work on all types implementing an informal interface (see information box below for more details).
This is similar to the philosophy employed in the C++ STL (Standard Template Library).

The following two lines make use of template programming to define two strings of type char, a text and a pattern.

.. includefrags:: demos/tutorial/a_first_example/basic_code_detailed.cpp
   :fragment: sequences

In order to store the similarities between the pattern and different text positions we additionally create a string storing integer values.

.. includefrags:: demos/tutorial/a_first_example/basic_code_detailed.cpp
   :fragment: score

Note that in contrast to the first two string definitions we do not know the values of the different positions in the string in advance.
In order to dynamically adjust the length of the new string to the text we can use the function :dox:`StringConcept#resize`.
The resize function is not a member function of the string class because SeqAn is not object oriented in the typical sence (we will see later how we adapt SeqAn to object oriented programming).
Therefore, instead of writing ``string.resize(newLength)`` we use ``resize(string, newLength)``.

.. includefrags:: demos/tutorial/a_first_example/basic_code_detailed.cpp
   :fragment: resize

.. note::

    Global function interfaces.

    SeqAn uses **global interfaces** for its data types/classes.
    Generally, you have to use ``function(variable)`` instead of ``variable.function()``.

    This has the advantage that we can extend the interface of a type outside of its definition.
    For example, we can provide a ``length()`` function for STL containers ``std::string<T>`` and ``std::vector<T>`` outside their class files.
    We can use such global functions to make one data type have the same interface as a second.
    This is called **adaption**.

    Additionally, we can use one function definition for several data types.
    For example, the alignment algorithms in SeqAn are written such that we can compute alignments using any :dox:`String` with any alphabet:
    There are more than 5 :dox:`String` variants in SeqAn and more than 8 built-in alphabets.
    Thus, one implementation can be used for more than 40 different data types!

After the string initializations it is now time for the similarity computation.
In this toy example we simply take the pattern and shift it over the text from left to right.
After each step, we check how many characters are equal between the corresponding substring of the text and the pattern.
We implement this using two loops; the outer one iterates over the given text and the inner loop over the given pattern:

.. includefrags:: demos/tutorial/a_first_example/basic_code_detailed.cpp
   :fragment: similarity

There are two things worth mentioning here: (1) SeqAn containers or strings start at position 0 and (2) you will notice that we use ``++variable`` instead of ``variable++`` wherever possible.
The reason is that ``++variable`` is slightly faster than its alternative, since the alternative needs to make a copy of itself before returning the result.

In the last step we simply print the result that we stored in the variable ``````score`` on screen.
This gives the similarity of the pattern to the string at each position.

.. includefrags:: demos/tutorial/a_first_example/basic_code_detailed.cpp
   :fragment: print

Refactoring
-----------

At this point, we have already created a working solution!
However, in order to make it easier to maintain and reuse parts of the code we need to export them into functions.
In this example the interesting piece of code is the similarity computation, which consists of an outer and inner loop.
We encapsulate the outer loop in function ``computeScore`` and the inner loop in function ``computeLocalScore`` as can be seen in the following code.

.. includefrags:: demos/tutorial/a_first_example/code_encapsulation.cpp
   :fragment: all

The function computeScore() now contains the fundamental part of the code and can be reused by other functions.
The input arguments are two strings.
One is the pattern itself and one is a substring of the text.
In order to obtain the substring we can use the function :dox:`SegmentableConcept#infix` implemented in SeqAn.
The function call ``infix(text, i, j)`` generates a substring equal to ``text[i ... j - 1]``, e.g. ``infix(text, 1, 5)`` equals "ello", where ``text`` is "Hello World".
To be more precise, infix() generates a :dox:`InfixSegment Infix` which can be used as a string, but is implemented using pointers such that no copying is necessary and running time and memory is saved.

Assignment 2
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Replace the code in your current file by the code above and encapsulate the print instructions.

   Hint
     The function head should look like this:

     .. container:: foldable

        .. includefrags:: demos/tutorial/a_first_example/solution_2.cpp
           :fragment: head

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/a_first_example/solution_2.cpp
           :fragment: all

        Your output should look like this:

        .. includefrags:: demos/tutorial/a_first_example/solution_2.cpp.stdout

The Role of References in SeqAn
-------------------------------

Let us now have a closer look at the signature of ``computeScore()``.

Both the text and the pattern are passed *by value*.
This means that both the text and the pattern are copied when the function is called, which consumes twice the memory.
This can become a real bottleneck since copying longer sequences is very memory and time consuming, think of the human genome, for example.

Instead of copying we could use **references**.
A reference in C++ is created using an ampersand sign (``&``) and creates an alias to the referenced value.
Basically, a reference is a pointer to an object which can be used just like the referenced object itself.
This means that when you change something in the reference you also change the original object it came from.
But there is a solution to circumvent this modification problem as well, namely the word **const**.
A ``const`` object cannot be modified.

.. important::

   If an object does not need to be modified make it an nonmodifiably object using the keyword ``const``.
   This makes it impossible to *unwillingly* change objects, which can be really hard to debug.
   Therefore it is recommended to use it as often as possible.

Therefore we change the signature of computeScore to:

.. includefrags:: demos/tutorial/a_first_example/solution_3.cpp
   :fragment: head

Reading from right to left the function expects two ``references`` to
``const objects`` of type ``String`` of ``char``.

Assignment 3
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Adjust your current code to be more memory and time efficient by using references in the function header.

   Hint
     The function head for ``computeLocalScore`` should look like this:

     .. container:: foldable

        .. includefrags:: demos/tutorial/a_first_example/solution_3.cpp
           :fragment: head_local

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/a_first_example/solution_3.cpp
           :fragment: all

        Your output should look like this:

        .. includefrags:: demos/tutorial/a_first_example/solution_3.cpp.stdout

Generic and Reusable Code
-------------------------

As mentioned earlier, there is another issue: the function computeScore only works for Strings having the alphabet ``char``.
If we wanted to use it for ``Dna`` or ``AminoAcid`` strings then we would have to reimplement it even though the only difference is the signature of the function.
All used functions inside ``computeScore`` can already handle the other datatypes.

The more appropriate solution is a generic design using templates, as often used in the SeqAn library.
Instead of specifying the input arguments to be references of strings of ``char`` s we could use references of template arguments as shown in the following lines:

.. includefrags:: demos/tutorial/a_first_example/solution_4_templateSubclassing.cpp
   :fragment: template

The first line above specifies that we create a template function with two template arguments ``TText`` and ``TPattern``.
At compile time the template arguments are then replace with the correct types.
If this line was missing the compiler would expect that there are types ``TText`` and ``TPattern`` with definitions.

Now the function signature is better in terms of memory consumption, time efficiency, and generality.

.. important::

   The SeqAn Style Guide

   The :ref:`SeqAn style guide <infra-contribute-style-cpp>` gives rules for formatting and structuring C++ code as well as naming conventions.
   Such rules make the code more consistent, easier to read, and also easier to use.

   #. **Naming Scheme**.
      Variable and function names are written in ``lowerCamelCase``, type names are written in ``UpperCamelCase``.
      Constants and enum values are written in ``UPPER_CASE``.
      Template variable names always start with 'T'.
   #. **Function Parameter Order**.
      The order is (1) output, (2) non-const input (e.g. file handles), (3) input, (4) tags.
      Output and non-const input can be modified, the rest is left untouched and either passed by copy or by const-reference (``const &``).
   #. **Global Functions**.
      With the exception of constructors and a few operators that have to be defined in-class, the interfaces in SeqAn use global functions.
   #. **No Exceptions**.
      The SeqAn interfaces do not throw any exceptions.

   While we are trying to make the interfaces consistent with our style guide, some functions have incorrect parameter order.
   This will change in the near future to be more in line with the style guide.

Assignment 4
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Generalize the ``computeLocalScore`` function in your file.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/a_first_example/solution_4.cpp

        Your output should look like this:

        .. includefrags:: demos/tutorial/a_first_example/solution_4.cpp.stdout

.. _oop-to-seqan:

From Object-Oriented Programming to SeqAn
-----------------------------------------

There is another huge advantage of using templates: we can specialize a function without touching the existing function.
In our working example it might be more appropriate to treat ``AminoAcid`` sequences differently.
As you probably know, there is a similarity relation on amino acids: Certain amino acids are more similar to each other, than others.
Therefore we want to score different kinds of mismatches differently.
In order to take this into consideration we simple write a ``computeLocalScore()`` function for ``AminoAcid`` strings.
In the future whenever 'computerScore' is called always the version above is used unless the second argument is of type String-AminoAcid.
Note that the second template argument was removed since we are using the specific type String-AminoAcid.

.. includefrags:: demos/tutorial/a_first_example/solution_4_templateSubclassing.cpp
   :fragment: subclassing

In order to score a mismatch we use the function ``score()`` from the SeqAn library.
Note that we use the :dox:`Blosum62` matrix as a similarity measure.
When looking into the documentation of :dox:`Score#score` you will notice that the score function requires a argument of type :dox:`Score`.
This object tells the function how to compare two letters and there are several types of scoring schemes available in SeqAn (of course, you can extend this with your own).
In addition, because they are so frequently used there are shortcuts as well.
For example :dox:`Blosum62` is really a **shortcut** for ``Score<int, ScoreMatrix<AminoAcid, Blosum62_> >``, which is obviously very helpful.
Other shortcuts are ``DnaString`` for ``String<Dna>`` (:ref:`sequence tutorial <tutorial-datastructures-sequences>`), ``CharString`` for ``String<char>``, ...

.. _template-subclassing:
.. tip::

   Template Subclassing

   The main idea of template subclassing is to exploit the C++ template matching mechanism.
   For example, in the following code, the function calls (1) and (3) will call the function ``myFunction()`` in variant (A) while the function call (2) will call variant (B).

   .. includefrags:: demos/tutorial/a_first_example/example_tempSubclassing.cpp

Assignment 5
""""""""""""

.. container:: assignment

   Type
     Application

   Objective
     Provide a generic print function which is used when the input type is not ``String<int>``.

   Hint
     Keep your current implementation and add a second function.
     Don't forget to make both template functions.
     Include ``<seqan/score.h>`` as well.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/a_first_example/solution_5.cpp

        Your output should look like this:

        .. includefrags:: demos/tutorial/a_first_example/solution_5.cpp.stdout

Tags in SeqAn
-------------

Sometimes you will see something like this:

.. includefrags:: demos/tutorial/a_first_example/base.cpp
      :fragment: seqan_tags

Having a closer look you will notice that there is a default constructor call (``MyersHirschberg()`` ) within a function call.
Using this mechanism one can specify which function to call at compile time.
The ``MyersHirschberg()`` `` is only a tag to determine which specialisation of the ``globalAligment`` function to call.

**If you want more information on tags then read on** otherwise you are now ready to explore SeqAn in more detail and continue with one of the other tutorials.

There is another use case of templates and function specialization.

This might be useful in a ``print()`` function, for example.
In some scenarios, we only want to print the position where the maximal similarity between pattern and text is found.
In other cases, we might want to print the similarities of all positions.
In SeqAn, we use **tag-based dispatching** to realize this.
Here, the type of the **tag** holds the specialization information.

.. tip::

   Tag-Based Dispatching

   You will often see **tags** in SeqAn code, e.g. ``Standard()``.
   These are parameters to functions that are passed as const-references.
   They are not passed for their values but for their type only.
   This way, we can select different specializations at **compile time** in a way that plays nicely together with metafunctions, template specializations, and an advanced technique called [[Tutorial/BasicTechniques| metaprogramming]].

   Consider the following example:

   .. includefrags:: demos/tutorial/a_first_example/example_tags.cpp

   The function call in line (3) will call ``myFunction()`` in the variant in line (1).
   The function call in line (4) will call ``myFunction()`` in the variant in line (2).

The code for the two different ``print()`` functions mentioned above could look like this:

.. includefrags:: demos/tutorial/a_first_example/example_tags_for_print.cpp


If we call ``print()`` with something different than ``MaxOnly`` then we print all the positions with their similarity, because the generic template function accepts anything as the template argument.
On the other hand, if we call print with ``MaxOnly`` only the positions with the maximum similarity as well as the maximal similarity will be shown.

Assignment 6
""""""""""""

.. container:: assignment

   Type
     Review

   Objective
     Provide a print function that prints pairs of positions and their score if the score is greater than 0.

   Hints
     SeqAn provides a data type :dox:`Pair`.

   Solution
     .. container:: foldable

        .. includefrags:: demos/tutorial/a_first_example/solution_6.cpp

        Your output should look like this:

        .. includefrags:: demos/tutorial/a_first_example/solution_6.cpp.stdout

Obviously this is only a toy example in which we could have named the two ``print()`` functions differently.
However, often this is not the case when the programs become more complex.
Because SeqAn is very generic we do not know the datatypes of template functions in advance.
This would pose a problem because the function call of function ``b()`` in function ``a()`` may depend on the data types of the template arguments of function ``a()``.

The Final Result
----------------

Don't worry if you have not fully understood the last section.
If you have -- perfect.
In any case the take home message is that you use data types for class specializations and if you see a line of code in which the default constructor is written in a function call this typical means that the data type is important to distinct between different function implementations.

Now you are ready to explore more of the SeqAn library.
There are several tutorials which will teach you how to use the different SeqAn data structures and algorithms.
Below you find the complete code for our example with the corresponding output.

.. includefrags:: demos/tutorial/a_first_example/final_result.cpp

Your output should look like this:

.. includefrags:: demos/tutorial/a_first_example/final_result.cpp.stdout
