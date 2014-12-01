How To: Implement custom subclasses of existing classes
-------------------------------------------------------

TOC

Overview
~~~~~~~~

SeqAn contains a number of abstract base classes, e.g.
seqan:Class.String and seqan:Class.Index, and different specializations
called template subclasses, e.g. seqan:"Spec.Alloc String" and
seqan:Spec.Index\_ESA. In the following we describe the abstract
interface of some of the base classes that has to be implemented by
every subclass.

Index
~~~~~

The set of functions and meta-functions that every index has to provide
can be seen [seqan:Class.Index here]. Some functions, e.g. indexRequire
and indexSupplied, are implemented by the base class and need not to be
reimplemented by subclasses. The following functions are abstract and
must be implemented by every index subclass:

+----------------+----------------------------------------------------+
| **Function**   | **Description**                                    |
+================+====================================================+
| clear          | Clears all index tables (called fibres)            |
+----------------+----------------------------------------------------+
| getFibre       | Returns a specific index table selected by a tag   |
+----------------+----------------------------------------------------+
| indexCreate    | Creates a specific index table.                    |
+----------------+----------------------------------------------------+

Typically, an index would provide more than only the construction
interface, e.g. searching a pattern, iterate it like a suffix tree.

Finder Interface
^^^^^^^^^^^^^^^^

Searching in SeqAn is provided by the seqan:Class.Finder interface. To
provide a search interface for an index, a Finder must be specialized
for this index. There already exists a generic specialization of the
Finder class for searching indices. For indices having a suffix array it
suffice to specialize the internal function \_findFirstIndex and to
define a default find algorithm via the meta-function DefaultFinder. See
the find interface of the q-gram index for an
[source:trunk/include/seqan/index/find\_index\_qgram.h example].

Suffix Tree Interface
^^^^^^^^^^^^^^^^^^^^^

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
