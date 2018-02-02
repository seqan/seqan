Yara - Yet another read aligner
===============================


Overview
--------

Yara is an *exact* tool for aligning DNA sequencing reads to reference genomes.

Main features
~~~~~~~~~~~~~

* Exhaustive enumeration of sub-*optimal* end-to-end alignments under the edit distance.
* Excellent speed, memory footprint and accuracy.
* Accurate mapping quality computation.
* Support for reference genomes consisiting of million of contigs.
* Direct output in SAM/BAM format.

Supported data
~~~~~~~~~~~~~~

Yara has been tested on DNA reads (i.e., Whole Genome, Exome, ChIP-seq, MeDIP-seq) produced by the following sequencing platforms:

* Illumina GA II, HiSeq and MiSeq (single-end and paired-end).
* Life Technologies Ion Torrent Proton and PGM.

Quality trimming is *necessary* for Ion Torrent reads and recommended for Illumina reads.

Unsupported data
~~~~~~~~~~~~~~~~

* RNA-seq reads spanning splicing sites.
* Long noisy reads (e.g., Pacific Biosciences RSII, Oxford Nanopore MinION).

Installation from sources
-------------------------

The following instructions assume Linux or OS X. For more information, including Windows instructions, refer to the `SeqAn getting started tutorial <http://trac.seqan.de/wiki/Tutorial/GettingStarted>`_.

Software requirements
~~~~~~~~~~~~~~~~~~~~~

**A modern C++11 compiler with OpenMP 3.0 extensions is required to build Yara. If unsure, use GNU G++ 4.9 or newer.**

* Git.
* CMake 3.2 or newer.
* G++ 4.9 or newer.

Download
~~~~~~~~

Yara sources are hosted on GitHub within the SeqAn library. Download the sources by executing:

::

  $ git clone https://github.com/seqan/seqan.git

Configuration
~~~~~~~~~~~~~

Create a build project by executing CMake as follows:

::

  $ mkdir yara-build
  $ cd yara-build
  $ cmake ../seqan -DSEQAN_BUILD_SYSTEM=APP:yara -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.9

Build
~~~~~

Invoke make as follows:

::

  $ make all

Installation
~~~~~~~~~~~~

Copy the binaries to a folder in your *PATH*, e.g.:

::

  # cp bin/yara* /usr/local/bin


Usage
-----

Yara consists of two executables:

* **yara_indexer** builds the index of a reference genome.
* **yara_mapper** maps DNA reads on the indexed reference genome.

This document explains only basic usage. To get complete usage descriptions, invoke each tool with -h or --help.

Indexer
~~~~~~~

Index a reference genome *REF.fasta.gz* by executing:

::

  $ yara_indexer REF.fasta.gz -o REF.index

**The indexer needs at least 25 times the space of the uncompressed reference genome**.
Be sure to dispose of that space inside the output folder.
The tool will take about one-two hours to index the human reference genome.
On success, the tool will create various files called *REF.index.**.

**The indexer does not work over GPFS and may have problems on other network filesystems**.

Mapper
~~~~~~

Single-end reads
^^^^^^^^^^^^^^^^

Map single-end DNA reads on the indexed reference genome by executing:

::

  $ yara_mapper REF.index READS.fastq.gz -o READS.bam

By default, the tool will report all co-optimal mapping locations per read within an error rate of 5%.
The results will be stored in a BAM file called *READS.bam*.

Paired-end reads
^^^^^^^^^^^^^^^^

Map paired-end reads by providing two DNA read files:

::

  $ yara_mapper REF.index READS_1.fastq.gz READS_2.fastq.gz -o READS.bam

Output format
^^^^^^^^^^^^^

Output files follow the `SAM/BAM format specification <http://samtools.github.io/hts-specs/SAMv1.pdf>`_.
In addition, Yara generates the following optional tags:

+-----+----------------------------------------------------+ 
| Tag | Meaning                                            | 
+=====+====================================================+ 
| NM  | Edit distance                                      |
+-----+----------------------------------------------------+ 
| X0  | Number of co-optimal mapping locations             |
+-----+----------------------------------------------------+ 
| X1  | Number of sub-optimal mapping locations            |
+-----+----------------------------------------------------+ 
| XA  | Alternative locations: (chr,begin,end,strand,NM;)* |
+-----+----------------------------------------------------+ 


Contact
-------

For questions or comments, feel free to contact: Enrico Siragusa <enrico.siragusa@fu-berlin.de>


References
----------

1. Siragusa, E. (2015). Approximate string matching for high-throughput sequencing. PhD Dissertation, Free University of Berlin.
2. Siragusa, E., Weese D., and Reinert, K. (2013). Fast and accurate read mapping with approximate seeds and multiple backtracking. Nucleic Acids Research, 2013, 1–8.
