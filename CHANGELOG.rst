SeqAn Changelog
---------------

This file summarizes the changes to the SeqAn library and apps.

Release 2.4.0
~~~~~~~~~~~~~

Library Features
^^^^^^^^^^^^^^^^

- Align
   - Generic parallelisation and vectorisation.
   - Support for SSE4, AVX2 and AVX512.
   - Speed-Ups of > 1000x compared to serial execution of long DNA alignments.
- Indexing
   - Parallel ``find()`` interface
   - Support for optimal search schemes (https://arxiv.org/abs/1711.02035)
- ReducedAminoAcid
   - several new Reduced Amino Acid alphabets are now available.
- VCF I/O
   - Now supports version 4.2 of the specification, i.e. columns with only eight fields.

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

- BAM I/O
   - Fix of jumpToRegion() functionality to start at the desired region instead of the index block.
   - Works correctly on big endian platforms now.
- Translation
   - Handle empty/too-short input correctly.
- Code Cleanup
   - various parts of the codebase have undergone cleanup and may now report deprecated functionality as being deprecated via compiler-functionality.

Platform Support
^^^^^^^^^^^^^^^^

- Compiler support:
   - Satisfies stricter warning levels of GCC-8, Clang-5 also with ``-std=c++17``.
   - Supports VS2017.
   - Intel Compiler suite 2016 **was dropped**.
   - Intel Compiler suites 2017 and 2018 are newly supported.
   - GCC-4.9 on MacOS **was dropped** (newer GCC available everywhere via MacPorts or Homebrew).
   - Clang-3.5 on FreeBSD **was dropped** (newer Clang available in base system and Ports).
- CPU architectures support:
   - Substantial fixes for big endian platforms
   - Now officially supported and passing integration tests:
      - ``i386, amd64/intel64, x32, ia64``
      - ``armel, armhf, arm64``
      - ``mips, mipsel, mips64el``
      - ``powerpc, ppc64, ppc64el``
      - ``s390x, alpha, m68k, sh4``
   - Officially **not** supported: ``sparc64``
   - Thanks to the `Debian Med team <https://www.debian.org/devel/debian-med/>`_ for their patches!
- Upstream packages:
   - SeqAn2 packages are finally coming to Fedora, thanks to @sagitter
   - Package updates in Debian, Ubuntu, MacPorts, Homebrew and FreeBSD expected shortly.

Release 2.3.2
~~~~~~~~~~~~~

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

- Argument parser
   - various fixes in the version checker
   - fix incompatibilities in CTD file creation with KNIME (introduced in 2.3.0)
- Build systems
   - reintroduce ``FindSeqAn.cmake`` for projects that rely on cmake's module mode
   - fix the pkgconfig file
- Platform related
   - improved compliance with warning levels of soon-to-be-released gcc7 and clang4
   - because of unresolved bugs we now recommend gcc5 as minimum gcc version when using static linking

Release 2.3.1
~~~~~~~~~~~~~

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

- Argument parser
    - bool option negative values
    - improve and fix version check

Release 2.3.0
~~~~~~~~~~~~~

Library Features
^^^^^^^^^^^^^^^^

- Argument Parser:
    - Adds version check support to the argument parser.
        - Check for new updates of a specific application.
        - Check for new versions of the library.
        - This option is opt-out by default but can be switched to opt-in or completely disabled via compiler flags and the SeqAn build system.
    - Altered Argument Parsers help page to display argument information.
    - Extended Argument types by bool, input_directory and output_directory.
    - Display file extensions that contain numbers.

- Sequence I/O:
    - New support for RNA structure files
        - Supported formats: Vienna (.dbv), Dot-Bracket-Notation (.dbn), Stockholm (.sth), Connect (.ct), Bpseq (.bpseq), Extended Bpseq (.ebpseq)
        - Input/output of whole files or of a single record/header
    - Added function isOpen() for formatted files.
    - Enabling assignment of format tags that differ from underlying format.
    - It is now possible to treat a BAM file as a (compressed) sequence file and read the sequences as if they were FastQ.

- Blast I/O:
    - Added support for handling the ``Q_ACC``, ``S_ACC``, ``S_ALLACC``, ``S_TAX_IDS`` fields
    - Added non standard fields ``LCA_ID`` and ``LCA_TAX_ID`` for lowest common ancestor information
    - Moved some redundant data from matches into record objects

- FM Index:
    - Added documentation for the bidirectional FM index
    - Reduced size of constant-time FM index

- Graphs:
    - Added new function getVertexAdjacencyVector()

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

  - Sequences:
      - Initialize empty CStyle Strings properly.
      - Fixed length function for const Dependent-StringSet

  - Graphs:
      - Reimplemented DFS in a non-recursive fashion to avoid stack overflow.
      - Multiple Sequence Alignment: Fix getAlignmentStatistics() on empty ``matches`` string.

  - Alignments:
      - Banded Chain Alignment: check for possible score overflow.

  - GFF / GTF:
      - Fixed I/O compatibility
          - Ignoring additional space
          - Allowing records to have multiple parents

  - BAM I/O:
      - Parsing the header for SO tags

  - VCF I/O:
      - Fixed reading of contig names in VCF header

  - Indices:
      - Enforce Container-Types for find()

App Updates
^^^^^^^^^^^

  - Gustaf:
      - Fixed name conflict (TANDEM)

Platform Support
^^^^^^^^^^^^^^^^

  - Compiler support:
      - SeqAn satisfies stricter warning levels of GCC7 and c++1z
  - New operating systems supported:
      - (Debian) GNU/kFreeBSD and GNU/Hurd
  - New CPU architectures supported:
      - ``arm`` and ``arm64``, ``mips`` and ``mips64``
      - ``powerpc``, ``powerpc64`` and ``sparc64``
      - and some others (all Debian platforms except ``sh4`` and ``armel``)
  - Thanks to the Debian Med team for their patches

Infrastructure Updates
^^^^^^^^^^^^^^^^^^^^^^

  - Added feature to selectively deactive the build of individual apps
  - Enforce using Python 2.x for documentation
  - Improvements to CMake and PkgConfig files

Release 2.2.0
~~~~~~~~~~~~~

Library Features
^^^^^^^^^^^^^^^^

- Indices:
    - FM index now has several options to reduce space consumption or improve running time
        - up to three level rank dictionaries
        - size of blocks on the lowest level (referred to as ``WORDS_PER_BLOCK``)
    - Bidirectional FM index with constant running time using EPR-dictionaries
    - Please see the `manual <seqan.readthedocs.io/en/master/Tutorial/DataStructures/Indices/FMIndex.html>`_ for more information

- Alignment:
    - Vectorized DP-Alignment algorithms using SSE3/AVX2. Allows for inter-parallel alignment computation in a many-vs-many or one-vs-many mode.
    - add a scoring matrix type that can be specified at runtime (e.g. BLOSSUM62, BLOSSUM50)

- Modifier:
    - ModifiedString ModPadding: Expand a string with padding symbols, without changing the source.

- Other:
    - Replace pthread implementation with STLs thread support library. Increases performance and fixes rare bugs in bam_io.

App Updates
^^^^^^^^^^^

- SAK (Swiss Army Knife):
    - fixed sequence filters.
- Yara:
    - verifying seeds
    - fixes CIGARs and secondary records.

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

- Alignments:
    - fixes MyersHirschberg implementation.
    - accept '=' operations in CIGAR string
- Split Alignment:
    - computes correct trace from split position.
    - allows flexible free-end gaps configuration.
- close Fasta file after FAI-Index is built.
- fixes Character to AminoAcid conversion.
- remove temporary files created during tests on Windows.

Infrastructure Updates
^^^^^^^^^^^^^^^^^^^^^^

- Build System:
    - The Intel Compiler is now fully supported on Linux and Windows, both 32bit and 64bit; it builds faster binaries and supports some functionality not available in MSVC.
    - On Windows there is now experimental support for Clang/C2, the Microsoft version of the clang compiler.
    - Please see the `manual <http://seqan.readthedocs.io/en/master/Infrastructure/Use/CMakeBuildDirs.html#visual-studio>`_ for more information on how to use these compilers.
    - support deb/rpm/exe/dmg packages and SSE4+POPCNT binaries

- Platforms:
    - full FreeBSD support
    - Ship UCRT, OPENMP and Intel DLLs for apps on windows
    - more apps available on Windows and some packaging fixes

Documentation Updates
^^^^^^^^^^^^^^^^^^^^^

- Api Docs:
    - Tree-View by Module

Release 2.1.1
~~~~~~~~~~~~~

Minor release including major improvements of the manual, several library bug-fixes and changes in the build system. All library modules are backward compatible
with 2.1.0. For a complete list of changes visit `GitHub <https://github.com/seqan/seqan/pulls?q=is%3Apr+is%3Amerged++milestone%3A%22Release+2.1.1%22+>`_.

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

- Tests:
    - delete automatically created temp directories in unit and app tests
    - demo tests: ``std::cout`` was not considered in tests

App Updates
^^^^^^^^^^^

- Yara:
    - fall back to single-end mapping when paired-end library length distribution is neither given nor estimable
    - fixed handling of reference metagenomes (references larger than 16k sequences)
    - enabled support for reference metagenomes by default (``-DYARA_LARGE_CONTIGS=ON``)
    - added option ``--sensitivity`` (low, high, full)
    - replaced option --output-secondary with ``--secondary-alignments`` (tag, record, omit)
    - renamed several options

Documentation Updates
^^^^^^^^^^^^^^^^^^^^^

- Manual:
    - major reworking of the manual
    - repaired links to API dox
    - hourly update of API dox for nightly builds

Infrastructure Updates
^^^^^^^^^^^^^^^^^^^^^^

- Build System:
    - more sensible execinfo detection
    - don't ship apps and the manual on library releases
    - introduce cmake ``-DSEQAN_OFFICIAL_PKGS=1`` to build upstream releases with static binaries
    - cache dependency detection on ``DEVELOP``
    - make it possible to do ``RELEASE_LIBRARY`` without dox

- Platforms:
    - basic BSD support
    - fixed warnings on Windows

- KNIME:
    - packaging - more flexibility when generating KNIME plugins of external apps


Release 2.1.0
~~~~~~~~~~~~~

Major release with many new features and applications.
Except where noted below, this release is compatible to previous 2.x releases.
For a complete list of changes visit `GitHub <https://github.com/seqan/seqan/pulls?q=is%3Apr+is%3Amerged++milestone%3A%22Release+2.1.0%22+>`_.

Library Updates and Selected Bugfixes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Apps:
    - Yara: fixed warnings, build errors and bugs, updated test files
    - Yara: new features (compute mapping qualities, estimate distribution of paired-end insert sizes)
    - Yara: follow SAM recommended practices for paired-end reads
    - T-Coffee: new feature ``deep coffee`` (aligning several hundred sequences)
    - Gustaf: introduced two phase breakpoint combination; updated readme and help messages
    - Removed old apps: Razers2, Flexbar and SeqCons in favor of newer releases

- Alignments:
    - added feature to count gaps to the left a of a position/iterator
    - disallow wrong use of scoring scheme for Hirschberg algorithm
    - extended AlignmentStats by number of gaps and length of the alignment
    - fixed evaluation of alignment
    - using gaps for integrateAlign and align_extend

- BLAST (new module):
    - E-Value statistics, including precomputed constants, bit-score and e-value calculation for alignments
    - support for reading and writing BLAST Tabular files (with and without comments)
    - support for writing BLAST Report files

- Indices:
    - added public function for trie and radix tree construction
    - Q-gram Index: allows sorting the hash-table according to the number of occurrences to reduce cache misses

- IO:
    - Tabix index: allowing range queries on chromosomal file formats such as VCF
    - Fai Index: optimized fasta index construction
    - BAM: added function to write tags from BamTagsDict to the tags field of a bam record
    - BAM: allowed BamTagsDict to take const CharStrings

- Misc:
    - fixed Iupac alphabet by replacing ``=`` by ``U``
    - added missing ``O`` character to amino acid alphabet
    - Argument Parser: a few new features such as help string for advanced options
    - removed random number engine and replaced it by the STL one
    - ZipIterator & ZipContainerView: iterating simultaneously over multiple containers
    - extended edges in graphs to store a reference to its source

- Modifier:
    - ModifiedString ModPos: iterating over a sequence in a predefined order
    - overload save() of ModifiedStrings for const strings
    - fixed Modified Iterators and ModView

- Journaled String Tree (new module):
    - reference compressed string set structure
    - for more details see the `publication <http://bioinformatics.oxfordjournals.org/content/30/24/3499.short>`_

- STL containers:
    - added a completely new adaptation to SeqAn interfaces that supports all STL containers, also ``std::array`` and ``std::forward_list``
    - greatly improved compatibility of SeqAn algorithms with STL containers so these can be used instead of SeqAn Strings

- Streams:
    - improved ZipStream

- Compatibility to previous versions
    - the random module was removed, please use the STL's random module instead
    - the ``StringSet<T, Dependent<Tight> >`` has been deprecated and will likely be removed for the next release
    - some SeqAn Macros have been deprecated since C++11 is now required, e.g. there is no ``SEQAN_AUTO_PTR_NAME``, only ``unique_ptr<>``
    - ``SEQAN_NAMESPACE_MAIN`` has been moved into the ``seqan`` namespace, so some of your Metafunction overrides may need to be adapted

Documentation Updates
^^^^^^^^^^^^^^^^^^^^^

- Dox:
    - added version selector in API dox


Infrastructure Updates
^^^^^^^^^^^^^^^^^^^^^^

- Build System:
    - Major improvements to build system resulting in cleanup and dropped dependencies
    - C++11 is now required and many datatypes now have move cosntructors and -assignment operators
    - added support for new compiler versions, but dropped support for older compilers
    - requirements are now GCC ≥ 4.9 or LLVM / Clang ≥ 3.5 (for Linux, Mac OSX, FreeBSD) and Visual C++ ≥ 14.0 / Visual Studio ≥ 2015 (for Windows)


Release 2.0.2
~~~~~~~~~~~~~

Minor release including several library bug-fixes as well as better documentation and infrastructure.
All library modules are backward compatible with 2.0.1.
For a complete list of changes visit `GitHub <https://github.com/seqan/seqan/pulls?q=is%3Apr+is%3Amerged++milestone%3A%22Release+2.0.2%22+>`_.

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

- Sequences:
    - fixed insert() for packed_string
    - fixed segfault bug for upac assignment in Dna5 StringSet
    - added insertValue(), insert() and replace() for StringSets
    - added empty() for std::list

- IO:
    - BAM I/O: adding spport for custom tags with floats
    - BAM I/O: BamTagsDict allows wrapping a const object
    - FastQ: fixed readRecord() for malformed fastq files (avoid skipping records)
    - FaiIndex: fixed readSequence/readRegion allocation

- Apps:
    - Gustaf: loading Fasta files with Iupac characters

Documentation Updates
^^^^^^^^^^^^^^^^^^^^^

- Dox:
    - fixed page redirection
    - minor bugs
    - code snippets in the documentation now undergo build tests and continuous integration to avoid outdated documentation

Infrastructure Updates
^^^^^^^^^^^^^^^^^^^^^^

- Platform Support:
    - FreeBSD support
    - updated prerequisites for GCC to >= 4.7 and Clang to >= 3.3
    - fixed warnings for gcc6
    - clang-3.7.x: deactivated openmp because of bug
    - fixed compiler-warnings in Visual Studio (/W2 produces no warnings anymore)
    - added support for Visual Studio 2014 and 2015

- Build System:
    - added pkg-config support
    - changed includes search priorities for CMake's FindSeqAn

- Continuous Integration:
    - added more platforms on TravisCI


Release 2.0.1
~~~~~~~~~~~~~

Minor release including several library bug-fixes as well as better documentation and infrastructure.
All library modules are backward compatible with 2.0.0.
For a complete list of changes visit `GitHub <https://github.com/seqan/seqan/pulls?q=is%3Apr+is%3Amerged++milestone%3A%22Release+2.0.1%22+>`_.

Library Bug Fixes
^^^^^^^^^^^^^^^^^

- Basic:
    - Added AminoAcid symbol "O"
    - Disabled global exception handler by default

- Sequence:
    - Added missing overloads for const Strings
    - Fixed and tested StringSet
    - Reworked STL containers adaption
    - Fixed several bugs in ModifiedStrings and ModifiedIterators

- Stream:
    - Worked around I/O with std::string
    - Supported multi-stream gzip files produced by Illumina Casava
    - Fixed BgzfStream tell()

- SeqIO:
    - Changed Raw file extension from .txt to .raw

- BAM I/O:
    - Fixed BIN computation
    - Fixed a bug in jumpToOrphans()
    - Fixed internal concurrency problems
    - Fixed readBamHeader() to clear the BamHeader
    - Added assertions to writeRecord()
    - Added BamIndex::save() to save .bai files

- Gff I/O:
    - Fixed parsing of comment lines

- FragmentStore:
    - Fixed loading Gtf/Gff3 files

- Index:
    - Fixed open() and save() for WT FMIndex
    - Added open() and save() for OpenAddressing QGramIndex

- Seeds:
    - Fixed a bug in sparse chaining
    - Fixed a bug in banded chain alignment

Documentation Updates
^^^^^^^^^^^^^^^^^^^^^

- Manual:
    - Fixed and improved several Tutorials and HowTos
    - Added version-aware links to the dox

- Dox:
    - Added @datarace entity
    - Fixed broken links in "See Also" section
    - Fixed a problem with close button in the side pane
    - Documented class VirtualStream

- Demos:
    - Restructured demos directory
    - Fixed several broken demos

Infrastructure Updates
^^^^^^^^^^^^^^^^^^^^^^

- Platform Support:
    - Added support for GCC 4.9 and Clang 3.7
    - Preliminary support for Clang 3.8 with OpenMP
    - Preliminary support for Visual Studio 2015
    - Preliminary support for FreeBSD 10.2

- Build System:
   - Tested all demos
   - Upgraged TravisCI builds to run on Docker
   - Fixed Java detection


Release 2.0.0
~~~~~~~~~~~~~

Major release with many new features and applications.
Note, the majority of the modules are backward compatible to the previous version.
Some modules, e.g. I/O-modules, have some adapted easier-to-use or unified interfaces.

Library Updates
^^^^^^^^^^^^^^^

- Faster and easier-to-use modules for basic and formatted file I/O:
    - ``stream``
    - ``seq_io``
    - ``bam_io``
    - ``vcf_io``
    - ``gff_io``
- Faster data structures:
    - FMIndex (up to 4X).
    - Packed Strings.
- New alignment modules:
    - X-Drop extension for alignments (``align_extend``)
    - Sequence-profile alignments (``align_profile``)
- New AminoAcid-Dna translation module (``translation``)
- The motif finding module (``find_module``) has been removed.

Infrastructure Updates
^^^^^^^^^^^^^^^^^^^^^^

- The repository has been migrated to GitHub (https://github.com/seqan/seqan).
- Continuous integration builds happen on TravisCI.
- The manual has been migrated to sphinx (http://seqan.readthedocs.io).
- The ``core`` and ``extras`` subfolders have been removed.

New Apps
^^^^^^^^

- ANISE and BASIL
    - Methods for the detection and assembly of inserted sequence in High-Throughput Sequencing Data.

- BS Tools
    - Bisulfite read mapping and SNP and methylation level calling.

- Fiona
    - A parallel and automatic strategy for read error correction.

- Gustaf
    - Generic mUlti-SpliT Alignment Finder.

- Mason 2
    - A read simulator.

- NGS ROI
    - Region of Interest Analysis for NGS Data.

- Samcat
    - Concatenate and convert SAM/BAM files (faster than samtools).

- Seqcons 2
    - Compute consensus from sequences sequences with and without approximate alignment information.

- Yara
    - Yet another read aligner (replaces Masai).


Release 1.4.2
~~~~~~~~~~~~~

Documentation-only release backward compatible with 1.4.1.


Release 1.4.1
~~~~~~~~~~~~~

This minor release should be backward compatible with 1.4. It contains small fixes and many demos for improving the API documentation. Some file format functionality has been added.

Highlights
^^^^^^^^^^

- Many new demos and improved API documentation throughout the library.
- New file format support and tutorials for this functionality: VCF I/O, BED I/O, and improvements to GFF and GTF I/O.

Selected Bug Fixes
^^^^^^^^^^^^^^^^^^

- ``gff_io.h`` does not contain corrupt includes any more
- Gapped X-drop seed extension now works with score matrices such as BLOSUM60.
- SAM writer code now writes ``255`` for invalid ``MAPQ`` and ``0`` for invalid/unapplicable ``TLEN`` instead of ``*``.
- Fix in Postorder ParentLinks VSTree Iterator.
- ``SEQAN_PATH_TO_ROOT()`` can now be used in demo programs.
- Removing duplicate definition of ``SEQAN_ENABLE_TESTING`` in build system.
- Write support for ``char *`` for ``BamTagsDict``.
- Fix in ``StringEnumerator``.
- Fix writing out of file extension when writing KNIME plugins.

Release 1.4
~~~~~~~~~~~

Highlights
^^^^^^^^^^

- New read mappers applications Masai and RazerS 3.
- Extended and more robust I/O functionality in ``stream``, ``seq_io``, ``bam_io``, and ``gff_io``.
- Module arg_parse creates improved command line help and supports workflow engine integration.
    - Also see https://github.com/genericworkflownodes
- Greatly improved alignment module with better performance and interfaces.
- Greatly improved build system, ``find_package(SeqAn)`` for your CMake build systems.

New Apps
^^^^^^^^

- ALF
    - Alignment free sequence comparison.

- Breakpoint Calculator
    - Breakpoint computation for genomic alignments.

- Masai
    - Fast index-based read mapper.

- RazerS 3
    - Fast filtration-based, parallel read mapper.

- SnpStore
    - SNP and small indel calling.

Major App Updates
^^^^^^^^^^^^^^^^^

- All applications now use the ArgumentParser and have better CLI help.

- Rabema
    - Rewritten from scratch, includes BAM support.
    - Greatly lowered memory requirements.

- SeqCons
    - Fixing input bugs, supports SAM I/O now.

- Stellar
    - Major update improving running time, including bug fixes, and
      allowing for various alphabet types.

- MicroRazerS
    - Adding support for SAM output.

Major Library Updates
^^^^^^^^^^^^^^^^^^^^^

- Modules ``seq_io``, ``bam_io``, ``gff_io`` with I/O functionality.
- FM Index in module ``index``.
- Rewritten ``align`` module with better performance, more consistent interfaces.
- Split alignment module ``align_split``.
- Metaprogramming: introducing ``EnableIf``, ``DisableIf``, ``EnableIf2``, and ``DisableIf2`` metafunctions
- Module ``alignment_free`` for alignment free sequence comparison.
- Module ``journaled_set`` for managing many similar sequences.
- Faster open addressing q-gram index.
- generic support for memory mapped files via FileMapping class
- Adding module ``parallel`` with atomic operations in C++98.
- Greatly improved FragmentStore documentation.
- Adding ``position()``, ``operator-()``, ``operator[]`` with proxy functionality and relation operators to journaled string iterator.
- Pigeonhole-based filter algorithm.
- Parallel repeat finding.
- Clang support, C++11 support

Major Library Bug Fixes
^^^^^^^^^^^^^^^^^^^^^^^

- Fixing repeat finding on Dna5Q.
- Fixing insert size computation in store_all.h
- Fixing memory initialization problem in ``appendValue()`` for Block String.
- Default constructor of Iter modified, such that data_container and data_position are initialized.
- Fixed error loading Fasta on Windows.
- Fixed wrong StringSet size types, allow to easily subclass Alloc strings
- Now supports SAM files with missing read sequences
- Fixing SeqAn code for C++11
- FragmentStore fixes.

Miscellaneous
^^^^^^^^^^^^^

- Experimental support added platforms for ICC and PGI compilers.
- Experimental support for CUDA.
- Build System
    - Large updates to build system.
    - Includes ``FindSeqAn.cmake`` for easily using SeqAn in your own CMake build system.
    - Packaging now based on CPack
- Xcode plugin for MacPorts LLVM/Clang in Xcode 3 and 4
- Improved code generator ``skel.py``.
- Many minor bug fixes
- Cleaned code base
- Added test cases (e.g. Stellar)
- Improved documentation and added examples (Mason, Rabema, RazerS, etc.)
- Improving coding style compliance of Array String implementation.
- Various tool improvements (e.g. RazerS 3)
- Performance improvements.
