.. sidebar:: ToC

    .. contents::

.. _infra-misc-write-nice-unix-programs:

Writing Nice Unix Programs
==========================

In bioinformatics, many programs are of "academic" quality, i.e. they are written to present a new method.
The implementation is often "only" used by other academics who, since they do not pay for the program, are willing to take some glitches for granted.

This page tries to help academic software developers in writing better Unix programs.
The points raised here come from our own experience with using academic software.
The focus is on C and C++ programming in a Unix (e.g. Linux, Mac Os X) environment.
The hints should also be applicable to other languages such as Java, and in some way also Windows.

Program Return Codes
--------------------

Rationale
^^^^^^^^^

The ``main()`` method of a program should be ``0`` if there were no errors and a value different from ``0`` otherwise.

Explanation & Reasoning
^^^^^^^^^^^^^^^^^^^^^^^

The ``main()`` function should return an integer indicating whether the program completed running successfully or not.
A value of ``0`` indicates that no error occurred, a value not equal to ``0`` indicates that an error occurred.
You might consider returning different values for different kinds of errors, e.g. ``2`` for I/O errors, ``3`` for logic errors and ``1`` as a catch-all for any all errors.

This makes it easy for a calling script/program to check whether the program executed successfully.

Example
^^^^^^^

The following program returns ``1`` if it receives any argument and ``0`` if this is not the case.

.. code-block:: cpp

    #include <cstdio>

    int main(int argc, char ** argv)
    {
      if (argc > 1) {
        fprintf(stderr, "I do not like arguments!\n");
        return 1;
      }

      return 0;  // Everything went smoothly!
    }

The following bash script calls programs and reacts to the return code.

.. code-block:: console

    #!/bin/sh

    # 1. Only success case.
    program arg1 arg2 && echo "success!"

    # 2. Only failure case.
    {|
    ! echo "failure"
    |}


    # 3. Handle success/failure case
    program arg1 arg2
    if [ "$?" ]; then
      echo "success"
    else
      echo "failure"
    fi

    # 4. Use case for separating cases
    # TODO

Links
^^^^^

*  `Error Level @ Wikipedia <http://en.wikipedia.org/wiki/Exit_status>`_

Assume Few Things About Paths
-----------------------------

Rationale
^^^^^^^^^

Do not assume anything on paths for (1) the program to reside in (2) temporary files or (3) the working directory.
Fixing the program install directory at configure time is OK.

Explanation
^^^^^^^^^^^

Most Unix programs are configured with a ``$PREFIX`` variable (e.g. setting ``--prefix=`` in the ``./configure`` script) and assume that all paths are relative to the given one.
For example, the Apache 2 web server reads its configuration from the director\ ``${PREFIX}/apache2``.
This is a reasonable assumption. Another reasonable assumption is that the current working directory is writeable.
However, temporary files should be stored in ``${TMPDIR}`` or ``/tmp`` (see the related section).

Non-reasonable assumptions are:

*  *The program is executed in the directory the binary resides in.*
   For example, program ``prog`` at path ``/path/to/prog`` should not assume that the working directory is ``/path/to`` when it is executed.
   Especially, do not assume that the directory the binary resides in is writeable.
   If your program is installed in ``/usr/bin``, this path is non-writeable for normal users on Unix.
*  A program *must* be in a given specific path fixed at *code writing time*.
   While it is reasonable for the user to give an install path at *configure-time*, the user should be able to install the program in any directory, including ``/opt``, his ``${HOME}`` directory or ``/some-weird-path/the/sys/admin/gave``.

Best practice is:

*  Use ``${TMPDIR}`` if available, fall back to ``/tmp``, for intermediate/temporary files.
*  Use reasonable defaults for result files, e.g. the path the input file resides in.
*  Allow the user to set an output directory.
*  Consider asking the user before overwriting result files when using defaults.

Example
^^^^^^^

Some programs create the result files in the current working directory.
This is not good practice, since the current working directory is *context* dependent.
While it is possible to use ``pushd`` and ``popd`` to use one directory per call to the program, it is much less error prone and more comfortable for the caller to specify the file on the comman dline.

Provide Good Defaults
---------------------

Rationale
^^^^^^^^^

Require as few parameters as possible, provide defaults or guess as many as possible.

Explanation
^^^^^^^^^^^

The more parameters are required in a program, the hard it gets to use.
For many parameters, default values can be given by the program's author.
Other parameters can be guessed depending on the input.

It should still be possible to override such value by command line parameters.

Example
^^^^^^^

The quality type of a FASTQ file can be guessed from the file contents very reliably by looking at the quality entries.
Nevertheless, the user should be able to override this by a command line parameter.

Positional vs. Named Arguments
------------------------------

TODO
^^^^

Provide all-in-one-go Variants of your program
----------------------------------------------

Rationale
^^^^^^^^^

While many program's steps might add to flexibility, a tool is easier to use if there is only one call to it.

Explanation
^^^^^^^^^^^

Some bioinformatics programs consist of many steps, e.g. (1) building an index (e.g. k-mer or suffix array) (2) perform a search, and (3) combine multiple search results to one.
While this might enable the flexible usage of the program it complicates its usage.
Please also provide a way to call your program that creates an output from the input files with one program call.

Example
^^^^^^^

For paired-end read mapping, the program *bwa* consists of multiple calls.

#. Call bwa to build an index on your genome.
#. Map the left-end reads, yielding a position file.
#. Map the right-end reads, yielding a position file.
#. Combine the two position files previously created.

While it is OK to first create an index file (this file can be used for many reads files), the last three steps could be combine into one umbrella command.
This would reduce the number of intermediate files and be much more comfortable for users.

Use ``stdout`` and ``stderr`` correctly
---------------------------------------

Rationale
^^^^^^^^^

The standard stream ``stdout`` is for the program's output while ``stderr`` is for logging and error messages.
It should be possible to redired ``stdout`` to an output file and ``stderr`` to a log file.
Use ``-`` as shortcuts for ``stdout`` and ``stderr``.

Explanation
^^^^^^^^^^^

In C/Unix programming ``stdout`` is for output to the user, ``stderr`` is for error messages and logging.
For example, when running daemons (e.g. web servers), the output to ``stderr`` ends up in log files.

If your program has only one input and one output file, it could accept the input from ``stdin`` by default and write to ``stderr``.
An example is the ``grep`` tool on Unix. You can specify different programs on the command line, however.

If you have program arguments for input and output files then you should use ``-`` for shortcuts to ``stdint`` and ``stderr``.
For example, a call to ``program --in-file - --out-file -`` would read from ``stdin`` and write to ``stdout``.

Example
^^^^^^^

*  When the program is called with wrong parameters, the return code should not be ``0`` and the help should be printed to ``stderr``.
*  When the program is called with a ``--help`` parameter, the return code should return ``0`` and the help should be printed to ``stdout``.

Allow specifying all file names through the command line
--------------------------------------------------------

TODO
^^^^

Do Not Require A Specific Working Directory
-------------------------------------------

Rationale
^^^^^^^^^

Do not require that the current working directory is in any relation to the directory containing the binary.

Explanation
^^^^^^^^^^^

Some programs must be called with ``./program``, e.g. the current working directory.
This makes it harder to use the program when
installed centrally and when multiple instances are called at the same time on the same file system.
This makes it harder to use in complex software pipelines.
Here, additional working directories and either symbolic links or copies of the program binary have to be created for each called instance.

Use ``$TMPDIR`` For Temporary Files, Fall Back to */tmp*
--------------------------------------------------------

Rationale
^^^^^^^^^

Use the value of the environment variable ``${TMPDIR}`` for temporary files.
If it is not set, use ``/tmp`` or ``/var/tmp``.

Explanation
^^^^^^^^^^^

On Unix, the canonical place to store temporary file is the value of the environment variable ``${TMPDIR}``.
If it is not set, then use ``/tmp`` or ``/var/tmp``.
``/tmp`` might be cleared during system reboots while ``/var/tmp`` is not cleared during system reboots but possibly rather depending on the file age.

Links
^^^^^

*  `TMPDIR @ Wikipedia <http://en.wikipedia.org/wiki/TMPDIR>`_

Misc Links
----------

*  `Heng Li's "Debugging Memory Problems" <http://lh3lh3.users.sourceforge.net/memdebug.shtml>`_ (Heng Li of BWA, samtools etc. fame)
