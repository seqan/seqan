.. _manual-main:

SeqAn Manual
============

Welcome to the manual pages for the SeqAn library!

SeqAn is a C++ template library for the analysis of biological sequences.
As such, it contains algorithms and data structures for

* string representation and their manipluation,
* online and indexed string search,
* efficient I/O of bioinformatics file formats,
* sequence alignments, and
* many, many more.


Requirements
------------

Well, as SeqAn is written in C++ it might not come as a surprise, that you should bring some basic knowledge about the C++ programming language.
We made quite an effort to not depend on other libraries, than the on-board tools already bring.
This means, to learn SeqAn it suffices in the beginning to only know about C++ and the `STL <http://en.cppreference.com/w/cpp>`_ which is the standard template library defined by the ISO C++ committee.
The rest will be discussed in the subsequent tutorials step by step.

Before we start, here is a strong advice!
If you are diving into C++ for the first time, because you are new to programming or switched from another programming language, then we recommend, you first stroll through the `C++ FAQs <https://isocpp.org/faq>`_ to acquaint yourself with C++.
There you can find many useful tips about C++ and get some further readings.
It also will introduce you to the paradigms that we used for designing this library.
In the Getting Started section we will introduce you to the design decisions of SeqAn and lay down some of the basic programming paradigms we follow to make this library so efficient.
If you never heard about these paradigms, dont't worry.
We will give you code examples, which you can try out on your own.
Never forgot, there is no better way to learn a new language or language feature, than to actually program with it.
So keep your fingers attached to the keyboard and let's start right away!

If you didn't install SeqAn yet, please follow the :ref:`User Guide <infra-use>` instructions to install SeqAn first.
After that you should continue with the tutorials.

.. hint::
    Please note, that although we try hard to provide a very comprehensive list of topics, it is not always possible to cover every angle of the library and its features.
    The tutorials are thought as a first place to start.
    If you are more experienced with SeqAn you can use the :ref:`API documentation <api-documentation>` in addition to search for specific functions or classes.

.. _manual-main-tutorials:

Tutorials
---------

The tutorial section is organized such that you can efficiently search for a specific topic you want to learn more about.
Each tutorial takes 30 to 60 minutes of your time for learning how to use SeqAn.
So buckle up and jump right into using SeqAn using our tutorials!

    :ref:`tutorial-getting-started`
        These articles are required for every one that is new to SeqAn.
        Take your time and study these documents thoroughly, as they describe the fundamental concepts and design decisions of the library.
        Everything else depends on these informations.

    :ref:`tutorial-datastructures`
        In the data structure tutorials we introduce you to the main data structures of this library and their usage.
        Beginners should start with the :ref:`Sequence tutorial <tutorial-datastructures-sequences-strings-and-segments>`, and then continue with the :ref:`Alingment tutorials <tutorial-datastructures-alignment>`.
        After that beginners should continue with the :ref:`Alignment Algorithm tutorials <tutorial-algorithms-alignment>`.

    :ref:`tutorial-algorithms`
        In this section we explain several different algorithms that are crucial for many bioinformatics applications.
        This includes pattern matching, dynamic programming algorithms for sequence alignments, seed extension and many more.
        Beginners that come from the tutorials about data structures should either continue with :ref:`Online Pattern Matching <tutorial-algorithms-pattern-matching-online>` or with the :ref:`DP Alignment Algorithms <tutorial-algorithms-alignment>`.

    :ref:`tutorial-io`
        On this page you will learn how to read/write and work with common bioinformatic file formats, such as FASTA, BAM, BED, VCF files, and more.
        Beginners should start with the :ref:`File I/O Overview <tutorial-io-input-output-overview>`.
        This tutorial introduces you to the basic I/O concepts and data structures.

    :ref:`tutorial-how-to`
        The how-to page is divided into Recipes and Use Cases.
        The former section gives you some useful hints about miscellaneous topics.
        The latter section describes how some use cases can be solved with SeqAn.
        Things presented here are for experienced SeqAn users.
        If you are a beginner, first have a look at the tutorials above.

    :ref:`tutorial-workflows-index`
        These tutorials teach you how to integrate your application into workflow engines like KNIME or Galaxy.

..    :ref:`tutorial-apps`
        A brief overview of some of our own applications developed with the SeqAn library.

Infrastructure
--------------

    :ref:`infra-use`
        These articles describe how to get SeqAn, how to use it in your application and explain things you need to consider when building.
        Everyone should read it.

    :ref:`infra-contribute`
        Anyone who wants to contribute code or documentation to SeqAn should read this.
        You will learn about the conventions and coding style.

    :ref:`infra-manage`
        These pages cover the structure of the SeqAn repository, the git workflow and explain release procedures.
        All SeqAn team members should read this; and also downstream package maintainers.

.. _api-documentation:

API Documentation
-----------------

The API documentation can be found :dox:`mainpage here`.

Partners
--------
   

    |intelLink|_
    
    .. |intelLink| image:: Intel-Logo-300x198.jpg
                        :scale: 70%
    .. _intelLink: http://www.intel.com/


    |deNBILink|_

    .. |deNBILink| image:: deNBI_Logo_rgb.jpg
                        :scale: 20%
    .. _deNBILink: https://www.denbi.de


.. toctree::
    :caption: Tutorials
    :name: tutorial
    :hidden:
    :maxdepth: 1
    :titlesonly:

    Tutorial/GettingStarted/index
    Tutorial/DataStructures/index
    Tutorial/Algorithms/index
    Tutorial/InputOutput/index
    Tutorial/HowTo/index
    Tutorial/Workflows/index
..    Tutorial/Applications/index

.. toctree::
    :caption: Infrastructure
    :name: infra
    :hidden:
    :maxdepth: 1
    :titlesonly:

    Infrastructure/Use/index
    Infrastructure/Contribute/index
    Infrastructure/Manage/index
    Infrastructure/Misc/index

.. toctree::
    :caption: Follow Us
    :name: follow
    :hidden:
    :maxdepth: 1
    :titlesonly:

    SeqAn Project <http://seqan.de>
    GitHub <http://github.com/seqan>
    Twitter <http://twitter.com/SeqAnLib>
    RSS Feeds <http://www.seqan.de/feed/>
    SeqAn Mailinglist <https://lists.fu-berlin.de/listinfo/seqan-dev#subscribe>

.. toctree::
    :hidden:
    :caption: Appendix
    :titlesonly:

    zreferences

.. Generated pages, should not appear

    Indices and tables
    ------------------

    * :ref:`genindex`
    * :ref:`search`
