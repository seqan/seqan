.. sidebar:: ToC

    .. contents::

.. _infra-manage-write-tutorials:

Writing Tutorials
=================


At the bottom, you can find a `Tutorial Template`_ for starting a new tutorial.

Conventions
-----------

Wiki Conventions
----------------

* Use only one line per sentence. This increases the readability of the sources.

Naming Conventions
------------------

* Use `headline capitalization <http://www.newsletterfillers.com/archives/grammar/capitalization_headline.htm>`_ for headlines.
* Use the tutorial's title as the file name (e.g. ``/wiki/Tutorial/NameOfYourTutorial.rst``).
* Assignments are numbered in the order they appear in a tutorial (e.g. ``Assignment 5``).
  Do not use a section relative numbering but an absolute one.
  If, e.g., the last assignment of section 1 was assignment 3, the first assignment of section 2 is assignment 4).
* Place the assignment's solutions inline.

Design & Layout Conventions
---------------------------

* Use back ticks (``````) to denote names of variables, functions, etc. (e.g. ````append```` results in ``append``).
* Use bold font (``**word**``) to denote key concepts.
* Use ``item`` and ``menu > sub menu > item`` to denote GUI entries and menu paths.
* Use the following markup to include source code
  ::

      .. includefrags:: demos/tutorial/alignment/alignment_msa.cpp
         :fragment: init

  where ``demos/tutorial/tutorial/alignment/alignment_msa.cpp`` gives the source code file in the repository and ``init`` the fragment to include in the tutorial.
* You should always build and test the tutorials code snippets before using them.

  .. code-block:: console

     manual # make html

* Use the following markup to format screen output:
  ::

      ::

          # Hello World!

* Use the following markup to inform about **important bugs** or other relevant issues.
  The content (and thereby the box itself) is always of **temporary** nature and should **only be used thriftily**.:

  ::

      .. warning::

         Warning goes here.

* Use the following markup to give **important information**.

  These boxes contain information that **should be kept in mind** since the described phenomenon is very likely to be encountered by the reader again and again when working with SeqAn.
  In contrast to the ``.. warning::``, this box type is of **permanent** nature and the given information are valid for a long time.

  ::

      .. important::

         Important information goes here...


  Use the following markup to give further / **optional information**.
  These are information that support the understanding but are too distinct to be put in a normal paragraph.:

  ::

      .. hint::

         Optional information goes here.

* Use the following markup to format assignments (for further details see `Assignments`_):

  ::

       .. container:: assignment

          The assignment goes here.

* Use ``:dox:`DocItem``` to create links to the SeqAn API dox documentation.

  .. important::

     Note that this will mereley generate the URLs that **dddoc** would create but does not perform any checking.
     Some examples:

     * :dox:`String`
       (``:dox:`String```)
     * :dox:`AllocString`
       (``:dox:`AllocString```)
     * :dox:`AllocString Alloc String`
       (``:dox:`AllocString Alloc String```)
     * :dox:`StringConcept`
       (``:dox:`StringConcept```)

Structure
---------

Meta Information
----------------

Place the directives for the side bar and the link target for the tutorial page directly before the tutorial title.

::

    .. sidebar:: ToC

       .. contents::


    .. _tutorial-datastructures-sequences:

    Sequences
    ---------


Based on the `Tutorial Template`_, provide information regarding:

learning objective
  Describe the learning objective in your own words.

difficulty
  Valid values: Very basic, Basic, Average, Advanced, Very advanced

duration
  In average how much time will a user spend on absolving this tutorial?
  If you expect more than 90 minutes please split your tutorial up into multiple ones.

prerequisites
  A list of absolved tutorials and other requirements you expect your reader to fulfill.

Introduction
------------

In the next paragraph introductory information are given that answer the following questions:

* What is this tutorial about?
* Why are the information important?
* What are the communicated information used for?
* What can the reader expect to know after having absolved the tutorial?

Section
-------

Introduction
^^^^^^^^^^^^

In each section's introduction part you answer the following questions:

* What is this section about?
* What are the central concepts in this section?
* What is your partial learning objective?

Explanations / Examples
^^^^^^^^^^^^^^^^^^^^^^^

The main part consists of the description of the topic.
This is the space where enough knowledge is transmitted to **enable the reader to solve all assignments**.
Further details are contained in the `Tutorial Template`_ and in the didactics section.

Try not to get lost in details.
If you have useful but still optional information to give use a ``.. note::`` directive.

Assignments
^^^^^^^^^^^

The assignments' purpose in general is to support the reader's understanding of the topic in question.
For this each assignment is of a special type (Review, Application and Transfer), has an objective, hints and a link to the complete solution.

Depending on the type of assignment the reader is guided through the assignment solving by providing him with partial solutions.

There must always be an assignments of type Review.
Assignments must always appear in an ascending order concerning their types and no "type gap" must occur.

Thus the only valid orders are:

* Review
* Review, application
* Review, application, transfer

The order Review, transfer is invalid since a "type gap" (application type missing) occurred.

All assignments must be accompanied by a solution.

Further Section
^^^^^^^^^^^^^^^

as many further sections as you like

Didactics
---------

Type
^^^^

As already mentioned in the assignment structure description each assignment is of one type.

These levels are

Review
  knowledge fortification (mainly through repetition, optionally with slight variations)

Application
  supervised problem solving (finely grained step-by-step assignment with at least one hint and the interim solution per step)

Transfer
  knowledge transfer (problem solving in a related problem domain / class)

Based on the chosen level you should design your assignment.

Duration
^^^^^^^^

The time needed to absolve a tutorial must not exceed 90 minutes.
Split your tutorial up (e.g. Tutorial I, Tutorial II) if you want to provide more information.


Language
^^^^^^^^

Make use of a simple language.
This is neither about academic decadence nor about increasing the learning barrier.
You are not forced to over-simplify your subject but still try to use a language that is also appropriate for those who don't fully meet the tutorials prerequisites.

Mental Model
^^^^^^^^^^^^

When your describe and explain your topic give as many examples as possible.
Try to adopt the reader's perspective and imagine - based on your target group and prerequisites - your reader's mental model.
The mental model can be described as an imagination of the interaction of central concepts.
Try to support the reader in developing a mental model that fits best to your topic.

Integration
-----------

* Add a link to your tutorial to ``Tutorial.rst`` and add a link to the ``.. toctree``.
* Above you stated the tutorials your tutorial has as prerequisites.
  Add the link in a way that all required tutorials are listed above your tutorial.

Tutorial Template
-----------------

::

    .. sidebar:: ToC

       .. contents::


    .. _tutorial-tutorial-template:

    Tutorial Template
    =================

    Learning Objective
      Describe the learning objective in your own words.
      **Example:**
      You will be able to write a tutorial that meets our quality standards.

    Difficulty
      [Very basic, Basic, Average, Advanced, Very advanced]
      **Example:**
      Basic

    Duration
      In average how much time will a user spend on absolving this tutorial?
      If you expect more than 90 minutes please **split your tutorial up** into multiple ones.
      **Example:**
      1 h

    Prerequisites
      A list of absolved tutorials and other requirements you expect your reader to fulfill.
      **Example:** :ref:`tutorial-getting-started-first-steps-in-seqan`, :ref:`tutorial-algorithms-pattern-matching`, English language

    This is the place where introductory need to be in given, e.g. "This page constitutes the template for all future SeqAn tutorials".

    Use this and optional further paragraphs to give the following information:

    * What is this tutorial about?
    * Why are the information important?
    * What are the communicated information used for?
    * What can the reader expect to know after having absolved the tutorial?

    .. warning::

       This is a warning message.

       Here you can inform users about important bugs or other relevant issues.

    Section
    -------

    Use this and optional further paragraphs to give the following information:

    * What is this section about?
    * What are the central concepts in this section?
    * What is your partial learning objective?

    When your describe and explain your topic give **as many examples as possible**.
    Try to adopt the reader's perspective and imagine - based on your target group and prerequisites - your **reader's mental model**.
    The mental model can be described as an imagination of the interaction of central concepts.
    Use a **simple language** and try to support the reader in developing a mental model that fits best to your topic.

    .. tip::

       What are tips for?

       An ``.. tip`` ist useful to give information that are **optional** and thus don't need to be read.
       Typical information are **further details** that support the understanding but are too distinct to be put in a normal paragraph.

       In this example you could tell the reader more about didactics and give him some useful links.

    .. important::

       What are important blocks for?

       These boxes contain information that **should be kept in mind** since the described phenomenon is very likely to be encountered by the reader again and again when working with SeqAn.

    Subsection
    ^^^^^^^^^^

    If you give code examples tell the reader what he can see and what is crucial to your snippet.
    Link all classes and other resources to the SeqAn documentation system by using ``:dox:Item` (e.g. :dox:`String`).
    In order to include code snippets use ``.. includefrags:: path``.

    .. includefrags:: demos/tutorial/alignments/alignment_banded.cpp
       :fragment: alignment

    If possible also include the generated output by given code in the console.
    Here is one example:

    .. code-block:: console

       0: ACAG
       1: AGCC
       2: CCAG
       3: GCAG
       4: TCAG

    Now that you gave an overview of important concepts of your topic let the user play with it!
    Formulate **small assignments** to allow the reader to fortify his newly acquainted knowledge.

    Assignment 1
    """"""""""""

    .. container:: assignment

       Type
         [Review, Application, Transfer]

         Note that your readers will be in different phases of learning. For the sake of simplicity we restrict ourselves to the following three levels:

         #. knowledge fortification (mainly through repetition, optionally with slight variations)
         #. supervised problem solving (finely grained step-by-step assignment with at least one hint and the interim solution per step)
         #. knowledge transfer (problem solving in a related problem domain / class)

         **Example:** Application

       Objective
         The objective of the assignment.
         **Example:**
         Output all symbols a given alphabet can have.
         The output should look like this: ...

       Hints
         ...

       Solution
         .. container:: foldable

            Foldable solution with description.

         This part of the assignment is to give partial solutions.
         A partial solution starts with a sentence of what this step is about and gives the lines of code that are needed to implement this step.

         Solution Step 1
           .. container:: foldable
             The given sequence are of alphabet...
             Therefore, you have to...

             .. includefrags:: demos/tutorial/alignments/alignment_banded.cpp
                :fragment: main

         Solution Step 2
           .. container:: foldable
             The given sequence are of alphabet...
             Therefore, you have to...

             .. includefrags:: demos/tutorial/alignments/alignment_banded.cpp
                :fragment: fragment
