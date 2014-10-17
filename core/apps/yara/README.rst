Yara - Yet another read aligner
===============================


Overview
--------

Yara is yet another read aligner.

Main features
`````````````

* Reporting of optimal end-to-end alignments under the edit distance
* Multi-core parallelization
* Low memory footprint
* Output in SAM or BAM format

Supported sequencing platforms
``````````````````````````````

* Illumina HiSeq and MiSeq
* Life Technologies Ion Torrent Proton and PGM

Reads quality control and trimming is necessary for Ion Torrent reads, yet 
recommended in general.

Supported sequencing protocols
``````````````````````````````

* DNA-seq
* Exome-seq
* ChIP-seq

Please, note that Yara does not support RNA-seq as it cannot map reads between 
intron-intron splicing sites.


Usage
-----

Yara consists of two executables:

* yara_indexer: builds an index of a reference genome;
* yara_mapper:  maps reads on an indexed reference genome.

This documents explain the most basic workflow. To get a complete usage 
description, call each program with -h or --help.

Indexer
```````

First index a reference genome, e.g. by executing:

::

  $ yara_indexer REF.fasta

The reference genome must be stored inside a DNA (multi-)Fasta file.
The indexer runs in about one-two hours on mammal reference genomes.

Important
'''''''''

The indexer needs a considerable amount of disk storage!

If the indexer fails to process the reference genome, please define a system
temporary folder of adequate capacity, either by setting the environment
variable TMPDIR, e.g.:

::

  $ TMPDIR=/big/folder/

or by passing the working temporary folder to the indexer, e.g.:

::

  $ yara_indexer --tmp-dir /big/folder/ REF.fasta


Mapper
``````

Now map DNA reads on the previously indexed reference genome by executing:

::

  $ yara_mapper REF READS.fastq

The mapper by default will report all co-optimal mapping locations per read 
within an error rate of 5%. Results will be stored in a SAM file called 
READS.sam.

To map paired-end reads, pass both paired-end read files to the mapper:

::

  $ yara_mapper REF READS_1.fastq READS_2.fastq

To increase the number of mapped reads, increase the error rate e.g. to 6%:

::

  $ yara_mapper --error-rate 6 REF READS.fastq


Website
-------

For more information and updates, visit: http://www.seqan.de/projects/yara.


Contact
-------

For questions or comments, feel free to contact:
  Enrico Siragusa <enrico.siragusa@fu-berlin.de>
