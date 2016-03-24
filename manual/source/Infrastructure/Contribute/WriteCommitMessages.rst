.. sidebar:: ToC

    .. contents::

.. _infra-contribute-git-commits:

Writing Commit Messages
=======================

Format
------

On every commit to our revision control system (currently SVN) please provide a commit message of the following form:

::

    [CLASS1,CLASS2,...] Short description

    Optional long description

*  The first line starts with an arbitrary number of tags in square brackets, e.g. ``[CLASS1]`` or ``[CLASS1,CLASS2]``.
   See below for a possible list of classes.
*  These tags are followed by a short description, try to keep the first line below 120 characters, 80 if possible.
*  You can add an optional long description after an empty line.


Possible Classes
----------------

| **NOP**
|    Only whitespace changes.
|    e.g. removed trailing whitespace, replaced tabs by spaces, changed indentation
| **DOC**
|    Changes in the user documentation.
|    This includes changes to the DDDoc documentation, README files etc.
| **COMMENT**
|    Changes in the source documentation.
|    These changes are not visible to the users.
|    This includes ``TODO(${name}):`` statements.

| **API**
|    Changes in the API.
|    These changes classically break backward compatibility.
|    e.g. renaming of function names, changing of function parameter order.
| **INTERNAL**
|    Changes in the implementation.
|    These changes do not influence the public the API.
|    e.g. renaming of variable names, simplification of code
| **FEATURE**
|    A user-visible feature.
|    e.g. extension of an interface, measurable performance improvement
|    *If the change is also API breaking the classes FEATURE* **and** *API must be used.*
| **FIX**
|    Bug removal.
|    If one or more bugs from the ticket tracker are removed then this should be written as ``[FIXED-#7,#35]`` where ``#7`` and ``#35`` are ticket numbers.
| **TEST**
|    Addition or changes of tests.
|    All code changes that are accompanied with tests must provide the original and the TEST class.
|    Don't consider this as a coercion but as a privilege to use both classes!
| **CLI**
|    Change to the command line interface of a program.
|    e.g. change to the arguments a program accepts, change to the messages printed to the user
|     *Output that is meant for debugging or detailed introspection is handled by the LOG class.*
| **LOG**
|    Change of output for developers or very advanced users.
|    This is the output that is meant for debugging or detailed introspection that is excluded from ``CLI``.
|    Such output is usually printed to ``stderr``.

Examples
--------

Example: API Changes
^^^^^^^^^^^^^^^^^^^^

API change with tests and detailed description of changes.

::

    [API,TEST] Large changes of align module's API.

    This is a large patch that mostly updates the align module:

    * The Anchor Gaps specialization is moved from the store module to the align module.
    * The Array Gaps class is rewritten.
    * Both Anchor and Array gaps have mostly the same API now, differences are documented.
    * Greatly unified the interface of the ``globalAlignment()`` and ``localAlignment()`` interface.
    * The ``LocalAlignmentFinder`` is now called ``LocalAlignmentEnumerator``.
    * Threw out unused DP algorithms (DP algorithm will be cleaned up in the future, see below).
    * Clipping of gaps works consistently in both Anchor and Array Gaps data structures.
    * Wrote a large number of tests for all changes.
    * Adjusted SeqAn library and apps to new API.

    All in all, this should make the module more usable, albeit breaking the interface in some cases.
    There will be a future change by Rene that strongly unifies the DP algorithms.
    This will not inflict another API change, though, and further simplify the align module.

Example: Bug Fixes
^^^^^^^^^^^^^^^^^^

A fix that solves two tickets:

::

    [FIX-#240,#356] Fixed iteration of ``ModifiedString``s.

    Quite involved fix that allows iteration of ``ModifiedString`` objects.

A fix that does not have a ticket:

::

    [FIX] Fixed reading of CIGAR string in module bam_io.

    There was a bug when reading the operation "F", which was translated to
    FLABBERGASTED.  Fixed this to the documented behavior.

Example: Internal Changes
^^^^^^^^^^^^^^^^^^^^^^^^^

An internal change, reordering of code without changing the public API.

::

    [INTERNAL] Reordering code in module sequence so no more generated forwards are needed.

An internal change might include tests and improved comments.

::

    [INTERNAL,TEST,COMMENTS] Greatly improved transmogrify module.

    Restructured the whole internal structure of the module, adding a large number of tests
    and improving the source-level documentation.  The user level documentation is still
    lacking and should be the target of a future change.

Example: Changes To Command Line Interface And Logging
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Changes to the command line interface:

::

    [CLI] Changed output of STELLAR such to unify scientific notation floats.

Changes to logging in an app:

::

    [LOG] Improved logging in RazerS 5.

    Much more detailed logging allows easier debugging.  Part of this should probably be
    commented out before the next stable release once the dust has settled and most
    bugs have been removed.

