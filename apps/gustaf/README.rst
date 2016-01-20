Gustaf
======

Generic mUlti-SpliT Alignment Finder
------------------------------------

Overview
--------

| Gustaf is a tool primarily designed for multi-split mapping of
  sequencing reads. Gustaf uses SeqAn’s exact local aligner Stellar to find partial read
  alignments. Using an exact dynamic programming approach, we refine the alignments
  around possible split positions to determine precise breakpoint locations at
  single- nucleotide level.
|  

| **Usage note**: Stellar is not a read mapper, and hence, Gustaf is not
  designed to replace any read mapper pipeline with SV detection on top. We
  recommend doing read mapping with Yara or your favourit read mapper and then first run
  Stellar and then Gustaf on the remaining unmappable reads.

Installation
------------

| Gustaf is distributed with SeqAn - The C++ Sequence Analysis Library
  (see http://www.seqan.de). Please check out the
  latest version of SeqAn and build a release version with:

::

    git clone https://github.com/seqan/seqan.git -b develop
    mkdir seqan-build
    cd seqan-build
    cmake ../seqan -DCMAKE_BUILD_TYPE=Release
    make gustaf

| Precompiled binaries (Linux 64-bit, Windows, Mac OS X) of Gustaf can
  be downloaded from the SeqAn project page: http://www.seqan.de/projects/gustaf/
  
Usage
-----

| Stellar is not a read mapper, and hence, Gustaf is not designed to replace any read mapper pipeline with SV detection on top. We  recommend doing read mapping with your favourit read mapper and then first run Stellar
  and then Gustaf on the remaining unmappable reads.

| To get a short usage description of Gustaf, you can execute
  ``gustaf -h`` or
 ``gustaf --help``.

::

    gustaf <GENOME FASTA FILE> <READ FASTA FILE> [Options]

| Here, Gustaf will first run Stellar internally on the given input files. Please have a look at the Stellar options and default values for match length, error rate, etc.

::

    gustaf <GENOME FASTA FILE> <READ FASTA FILE> -m <GFF MATCH FILE> [Options]

| Using the ``--matchfile (-m)`` option, Gustaf skips running Stellar and directly uses the matches, preferable precalculated using Stellar,  from the given GFF file. We recommend using this option for larger data  sets.

::

    gustaf <GENOME FASTA FILE> <MATES1 FASTA FILE> <MATES2 FASTA FILE> -m <GFF MATCH FILE> [Options]

| Given two read pair input files, Gustaf is run in paired-end mode. Note that option ``-m`` is mandatory in this case. Please also use  Gustaf’s helper app ``gustaf_mate_joining`` (see below) to join both read files  before calling Stellar.
|  
| Remarks: There are different versions paired-end (or mate pair) data can come in. Gustaf always expect two files, and in default expects the second mates (MATES2) to be reverse complemented in relation to the first mates (MATES1), see option ``-rc`` to turn this behaviour off.

Non-optional arguments
^^^^^^^^^^^^^^^^^^^^^^

| Gustaf always expects both a database and a query file in Fasta format.
| Important: Following conventions, all Ids (line starting with ‘>’) have to be unique already within the first part until the first whitespace! For example:

::

    >Gustaf|group=6|reads=9 readId=1
    >Gustaf|group=6|reads=9 readId=2

are not unique.

| Without any additional parameters, Gustaf would call Stellar and then chain those matches of each read, that have either an overlap of at least 0.5 (50% of each match length), a gap in the read between the matches of at  most 10 bp, and that miss at most 15 bp at the end or beginning of the read. The program calls Stellar with default options (see Stellar Readme).
|  
| Two output files will be generated: breakpoint.gff includes all structural variant breakpoints in GFF format, breakpoint.vcf in VCF format.
|  
| The default behaviour can be modified by adding options to the command line. See ``gustaf --help`` and ``gustaf_mate_joining --help``

Output Formats
^^^^^^^^^^^^^^

| Gustaf currently supports the GFF and the VCF output format for reporting breakpoints. Using option ``-do``, Gustaf creates a Dot file of each(!) read that can be converted to png and show a graph with all alignments used for the read.
|  
| The General Feature Format is specified by the Sanger Institute as a tab- delimited text format http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
|  
| See  http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41 for information about the VCF file format specifications.

Example
-------

| See the *example* folder for example runs of Gustaf. There is a genome file *adeno.fa* and a modified genome file *adeno\_modified.fa*, from the latter we created the read file *adeno\_modified\_reads.fa* for single-end mode and the files *adeno\_modified\_reads\_mates1* and *adeno\_modified\_reads\_mates2* for paired-end mode.

Single-end Mode
^^^^^^^^^^^^^^^

| There is a pre-calculated Stellar output file *stellar.gff* for the single-end read file computed by calling

::

        stellar adeno.fa adeno_modified_reads.fa -l 30 -o stellar.gff

The default calls for Gustaf would then be

::

        gustaf adeno.fa adeno_modified_reads.fa -st 1 -l 30 -gff st1_l30.gff \
          -vcf st1_l30.vcf

or

::

        gustaf adeno.fa adeno_modified_reads.fa -st 1 -m stellar.gff \
          -gff st1_l30_m.gff -vcf st1_l30.vcf

| Both calls produce an output file containing the same breakpoints. In the first run, Gustaf internally calls Stellar with parameter ``-l 30``. In the second run, Gustaf used the pre-calculated file with Stellar matches (use this option for larger datasets where you want to run Stellar separately and reuse the Stellar output for multiple Gustaf runs).

Paired-end Mode
^^^^^^^^^^^^^^^

| In paired-end mode, we join both read pair files before calling Stellar. This can be done using the app ``gustaf_mate_joining`` by calling

::

        gustaf_mate_joining adeno_modified_reads_mates1.fa adeno_modified_reads_mates2.fa \
            -rc -o adeno_modified_reads_joinedMates.fa

| There is a pre-calculated Stellar output file *stellar\_joinedMates.gff* for the paired-end read file *adeno\_modified\_reads\_joinedMates.fa* computed by calling

::

       stellar adeno.fa adeno.fa adeno_modified_reads_joinedMates.fa -l 30 \
            -o stellar_joinedMates_l30.gff

The Gustaf call would then look like this

::

        gustaf adeno.fa adeno_modified_reads_mates1.fa adeno_modified_reads_mates2.fa \
            -m stellar_joinedMates_l30.gff -st 1 -ll 1000 -le 30 -rc \
            -gff gustaf_adeno_pairedend_ll1000le30.gff \
            -vcf gustaf_adeno_pairedend_ll1000le30.vcf

| Note the ``-rc`` parameter. In this simulated data, mate1 and mate2 have the same orientation, so we prevent the second file from beeing reverse complemented.

Joining Paired-end Reads With Gustaf\_mate\_joining
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The Gustaf directory includes another app called gustaf\_mate\_joining that can be build using the ``make gustaf_mate_joining`` command. Gustaf\_mate\_joining is a small app that helps prepare paired-end data for usage with Gustaf to  incorporate paired-end information. This simple program takes as input two mate pair or paired-end files and outputs a file where both mate sequences have been joined together. The FASTA file with joined mates is an required input file for the paired-end mode of Gustaf. The tool assumes the mates in the second file to be reverse complemented compared to the first file. This behaviour can be turned off using the command line argument ``-rc``.
|
| Given only one input file and two output files, the program will split the reads from the input files at half length, and write the first half of each sequence as mates1 into the first output file and the reversed complemented second half of each sequence as mates2 into the second output file. Reverse complementing the sequences can again be turned off using ``-rc``.

|  To prepare the joined mate file for the paired-end example above, call

::

        gustaf_mate_joining adeno_modified_reads_mates1.fa \
            adeno_modified_reads_mates2.fa -rc -o adeno_modified_reads_joinedMates.fa

| The mates in this small example are both from the same strand, so we avoid reverse complementing the second input file by using ``-rc``.


Contact
-------

| For questions or comments, contact:
|  Kathrin Trappe kathrin.trappe@fu-berlin.de
