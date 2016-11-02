SAM Cat
========

Concatenate and convert SAM/BAM files quickly.

(c) Copyright in 2015 by David Weese.

Usage
-----

This tool reads one or many SAM/BAM files and outputs their concatenation as a SAM/BAM file. 
The headers of more than one input file will be merged.
All formats will be detected from file extensions or input streams.

```
samcat one.bam two.sam three.bam -o concat.bam
```

Input and output formats can be arbitrary, i.e. to convert a SAM into a BAM file use:

```
samcat input.sam -o output.bam
```

If the output file name is omitted the result is written to stdout in SAM format, use ```-b``` to select BAM.
If ```-``` is used as input file name, the input is read from stdin.

```
cat input.sam | samcat -b - > output.bam
```

Additional Information
----------------------

As the BAM format is more compact than SAM and faster to parse the best performance will be reached with BAM files.

Contact
-------

For questions or comments, contact David Weese <dave.weese@gmail.com>.
