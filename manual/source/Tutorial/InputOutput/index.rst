.. We create this roles for putting the "Introduction: etc. headings
    on this page without them displaying in the ToC.  This would break
    rendering the ToC correctly on readthedocs style.  The rubric
    headings are formatted using CSS.

.. role:: rubric-heading1
    :class: rubric-heading1
.. role:: rubric-heading2
    :class: rubric-heading2

.. _tutorial-io:

Input/Output
============

.. toctree::
    :hidden:
    :titlesonly:

    FileIOOverview
    SequenceIO
    IndexedFastaIO
    BlastIO
    SamAndBamIO
    VcfIO
    BedIO
    GffAndGtfIO


SeqAn supports many standard file formats used in bioinformatics applications.
The :ref:`first tutorial <tutorial-io-input-output-overview>` shows you the basic concepts and data structures for handling file I/O.
So this would be a good starting point to learn more about it.
Otherwise, feel free to look into the tutorials for your desired file format and start reading and writing your files.