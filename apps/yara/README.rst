Yara - Yet another read aligner
===============================


Overview
--------

Yara is an *exact* tool for aligning DNA sequencing reads to reference genomes.

Main features
~~~~~~~~~~~~~

* Exhaustive enumeration of sub-*optimal* end-to-end alignments under the edit distance;
* Alignment of single-end, paired-end and mate-pair reads;
* Speed of 10-50 Gbp/h on a desktop computer;
* Fine-grained multi-threading;
* Low memory footprint via a generalized FM-index;
* Direct output in SAM or BAM format.

Sequencing data
~~~~~~~~~~~~~~~

Yara has been tested with DNA reads (WGS, Exome and ChIP-seq) produced by the following sequencing platforms:

* Illumina GA II, HiSeq and MiSeq;
* Life Technologies Ion Torrent Proton and PGM.

Quality trimming is *necessary* for Ion Torrent reads, yet strongly recommended in general.
Note that Yara cannot map RNA-seq reads spanning splicing sites.


Installation
------------

For installation instructions and updates, visit: http://www.seqan.de/projects/yara.


Usage
-----

Yara consists of two executables:

* **yara_indexer** builds the index of a reference genome;
* **yara_mapper** maps DNA reads on the indexed reference genome.

This documents explain the most basic workflow.
To get a complete usage description, invoke each tool with -h or --help.

Indexer
~~~~~~~

Index a reference genome *REF.fasta.gz* by executing:

::

  $ yara_indexer REF.fasta.gz -o REF.index

**The indexer needs at least 25 times the space of the uncompressed reference genome**.
Be sure to dispose of that space inside the output folder.
The tool will take about one-two hours to index the human reference genome.
On success, the tool will create various files called *REF.index.**.

Mapper
~~~~~~

Single-end reads
^^^^^^^^^^^^^^^^

Map single-end DNA reads on the indexed reference genome by executing:

::

  $ yara_mapper REF.index READS.fastq.gz -o READS.bam --threads 16

The mapper will report all co-optimal mapping locations per read within an error rate of 5%.
he results will be stored in a BAM file called *READS.bam*.
The tool will use 16 working threads.

Paired-end / mate-pair reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Map paired-end or mate-pair reads by providing two DNA read files:

::

  $ yara_mapper REF.index READS_1.fastq.gz READS_2.fastq.gz -o READS.bam \
                --library-length 300 --library-error 200

Be sure to provide the expected insert distribution using *--library-length* and *--library-error*,
as well as the correct orientation using *--library-orientation*.

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
| XT  | Type: Unique/Repeat                                |
+-----+----------------------------------------------------+ 
| XA  | Alternative locations: (chr,begin,end,strand,NM;)* |
+-----+----------------------------------------------------+ 


Contact
-------

For questions or comments, feel free to contact: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
