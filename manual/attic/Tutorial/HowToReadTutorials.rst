How To Read Tutorials
---------------------

TOC()

This tutorial explains how to read tutorials and the used typographic
conventions.

General Tutorial Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~

A tutorial starts with a section that gives some meta information about
the tutorial:

| ``Learning Objective ::``
| `` A short description of what this tutorial is about.``
| ``Difficulty ::``
| `` How advanced is the presented material?``
| `` One of Very basic, Basic, Average, Advanced, Very advanced.``
| ``Duration ::``
| `` An estimate of how long it takes to complete the tutorial.``
| ``Prerequisites ::``
| `` A list of tutorials you should have read before and other knowledge.``

Then, it gives a textual overview of the presented material. The
material is then broken into sections. Mostly, it makes sense to read
them sequentially but you might also try to jump to what interests you
most and then backtrack in case you are missing some information.

The last section in all tutorials is **Further Steps** that lists
suggestions for how to continue after completing the tutorial. At the
bottom of the page, there also is a link to send us suggestions if you
have any.

Typographic Conventions
~~~~~~~~~~~~~~~~~~~~~~~

-  GUI elements are shown in *italic*, e.g. a menu point is described by
   *Start -> All Programs -> Accessories -> Command Prompt*.
-  File names are written in ``typewriter text``.
-  Where appropriate we use
   `links <http://en.wikipedia.org/wiki/Hyperlink>`__, e.g. to the
   documentation: seqan:Class.ArgumentParser.

Source code is usually written with highlighting in blocks as follows
but you might also see small snippets such as
``int main(int, char const **)`` using typewriter font.

::

    #cpp
    #include <iostream>

    int main(int argc, char const ** argv)
    {
        std::cout << "argc == " << argc << '\n';

        return 0;
    }

Command Line Output
~~~~~~~~~~~~~~~~~~~

Sometimes, we show some command line output, for example

::

    #ShellBox
    # ls
    build      CMakeLists.txt  docs   extras   misc    sandbox
    CHANGELOG  core        docs2  LICENSE  README  util

Lines starting with a hash show the command line where the user types in
the command.

Generally, we do not give full paths to the binary. For example, if we
just compiled the example application ``demos/graph_algo_bfs.cpp``
and call it, the shell would use ``graph_algo_bfs`` instead of
``./demos/graph_algo_bfs`` in the call.

::

    #ShellBox
    # graph_algo_bfs
    Adjacency list:
    0 -> 4,1,
    1 -> 5,0,
    2 -> 3,6,5,
    3 -> 7,6,2,
    ...

In the Getting Started tutorial, the command line usage is more complex
than in the rest of the tutorials (which focus on C++, not on work on
the command line). Here, we will give additional information such as the
current working directory.

For example, a Mac/Linux shell output could look as follows. The name of
the current directory is given before the ``#``:

::

    #ShellBox
    username # cd Development
    Development # cd /var/log
    log # ls

For Windows command line prompts, we will give the full to the current
working directory.

::

    #ShellBox
    C:\Users\username> cd C:\
    C:\> mkdir Development
    C:\> cd Development
    C:\Development> dir

Boxes
~~~~~

Throughout the tutorial, you will find the following boxes.

::

    #WarningBox
    '''Warning:''' Warning Box Title

    This is a box with a warning.
    It is used rarely, mainly to indicate that there exists a known issue at the place the box is located at.

::

    #ImportantBox
    This is a box with a very important message.

::

    #InfoBox
    '''Information:''' Information Box Title

    This box contains additional, in-depth information that you might want to refer back to in the future.

::

    #AssignmentBox
    '''Assignment X:''' This is an assignment.

    Solve the tasks in the assignments to intensify your learning of SeqAn.

     Type ::
      The type of the task, one of Reproduction, Application, Transfer.
      ''Reproduction'' assignments only require you to do small changes or copy and adjust small pieces of code.
      ''Application'' assignments require you to apply the material just explained in a similar situation.
      ''Transfer'' assignments require some kind of problem solving skill to transfer what you just learned to a different application.
     Objective ::
      A short description of the task's objective.
     Solution ::
      A link to the solution.
      The solution might also be hidden.
      Click ''more...'' to see an example.
     Hints ::
      Some tasks contain additional hints.

    <pre>#FoldOut
    ----
    You found the hidden solution!

.. raw:: html

   </pre>

Further Steps
~~~~~~~~~~~~~

-  Go to the `Tutorial Table Of Contents <Tutorial>`__

*' Submit a comment*'

If you found a mistake or have suggestions about an improvement of this
page press:
[/newticket?component=Documentation&description=Tutorial+Enhancement+for+page+http://trac.seqan.de/wiki/Tutorial/HowToReadTutorials&type=enhancement
submit your comment].

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
