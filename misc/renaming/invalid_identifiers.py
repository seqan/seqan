#!/usr/bin/env python
import sys
import re

from helpers import *

PROGRAM_USAGE = """
SeqAn invalid identifiers detection script.

USAGE: invalid_identifiers.py BASE_PATH

BASE_PATH is the root path of all the folders to be searched.
This script generates a list of invalid identifiers found in the code base,
paired with their suggested replacement string in the format ``"old: new"``,
one identifier per line.
The result is written to the standard output.
""".strip()

INVALID_IDENTIFIER = re.compile(r'\b_[A-Z_]\w*\b')
REPLACEMENT_ID = re.compile(r'\b(__?)(\w*)\b')
# The following IDs are exempted from replacement since they are either defined
# by some compiler (-specific library) or are solely used within a string.
VALID_IDENTIFIERS = map(
        lambda rx: re.compile(rx),
        [ '___+',
          '^__$',
          '_N',
          '_L',
          '_H',
          '__u?int64(_t)?',
          '_A123456',
          '__OPTIMIZE__',
          '__gnu_cxx',
          '_Resize_String', # will be done manually
          '_Fill_String',   #
          '_Transcript_',
          '_Confidence_99',
          '_PARSER_H',
	  '_POSIX_TIMERS',
          '__GNUC_MINOR__',
          '_S_IREAD',
          '_S_IWRITE',
          '_O_BINARY',
          '_O_CREAT',
          '_O_RDONLY',
          '_O_RDWR',
          '_O_TEMPORARY',
          '_O_TRUNC',
          '_O_WRONLY',
          '_KMER_H',
          '_MSC_EXTENSIONS',
          '_GLIBCXX_PARALLEL',
          '_FILE_OFFSET_BITS',
          '_POSIX_SYNCHRONIZED_IO',
          '__cplusplus',
          '__(force)?inline(__)?',
          '__alignof(__)?',
          '__attribute__',
          '__GLOBAL__',
          '_DELETIONS____',
          '_INSERTS______',
          '_REPLACEMENTS_',
          '__int128',
          '__SSE2__',
          '__m128i',
          '__VA_ARGS__',
          '__FILE__',
          '__LINE__',
          '__GET_OPT_H__',
          '_OPENMP',
          '__SINIX__',
          '__sgi',
          '__BEOS__',
          '__aix__',
          '__ICC',
          '__WATCOMC__',
          '__ADSPBLACKFIN__',
          '_BEOS',
          '__SUNPRO_CC?',
          '__tru64',
          '__FreeBSD__',
          '__ultrix',
          '__OPENBSD',
          '_MPRAS',
          '_HAIKU',
          '_SGI_COMPILER_VERSION',
          '_POSIX_C_SOURCE',
          '_XOPEN_SOURCE',
          '__OpenBSD__',
          '__AIX__',
          '__ADSP21000__',
          '__HAIKU__',
          '__riscos__',
          '__hpux',
          '__HP_aCC',
          '__riscos',
          '__hpua',
          '__GNUC__',
          '_ULTRIX',
          '_SCO_SV',
          '__DECCXX',
          '_XENIX',
          '__sgi__',
          '_WIN32',
          '__PGI',
          '__QNX__',
          '__APPLE__',
          '__AIX',
          '_SGI',
          '_AIX',
          '__XENIX__',
          '__INTEL_COMPILER',
          '__osf',
          '__linux__',
          '__sinix__',
          '__bsdos__',
          '__ADSPTS__',
          '__sun',
          '__sinix',
          '__NetBSD',
          '__FreeBSD',
          '__osf__',
          '__ultrix__',
          '__COMPILER_VER__',
          '__COMO__',
          '__linux',
          '__UNIX_SV__',
          '__HAIKU',
          '__WIN32__',
          '__NetBSD__',
          '__CYGWIN__',
          '_COMPILER_VERSION',
          '__BORLANDC__',
          '__TRU64__',
          '__MINGW32__',
          '__aix',
          '__BeOS',
          '__QNXNTO__',
          '__hpux__',
          '__IBMCPP__',
          '__IAR_SYSTEMS_ICC__',
          '__18CXX',
          '__HP_cc',
          '__SUNPRO_C',
          '__DECC',
          '__IBMC__',
          '_MSC_VER' ])

def valid(id):
    """
    Returns whether the given ``id`` is in fact valid and shouldn't be replaced.
    """
    return any(VALID_ID.match(id) for VALID_ID in VALID_IDENTIFIERS)

def find_all(file):
    """
    Returns all invalid identifiers found in a given ``file``.
    """
    f = open(file, 'r')
    result = []
    for line in f:
        matches = INVALID_IDENTIFIER.findall(line)
        invalids = [match for match in matches if not valid(match)]
        result += invalids

    return result


def replacement(orig):
    """
    Returns the replacement string for a given invalid identifier.
    """
    return REPLACEMENT_ID.sub(r'\2\1', orig)


def generate_replacements(ids):
    """
    Generates a dictionary of replacement strings for a list of invalid
    identifiers.
    """
    return dict([(original, replacement(original)) for original in ids])


def main():
    if len(sys.argv) != 2:
        print >>sys.stderr, 'ERROR: Invalid number of arguments.'
        print >>sys.stderr, PROGRAM_USAGE
        return 1

    results = {}
    project_path = sys.argv[1]

    for file in all_files(project_path):
        results[file] = set(find_all(file))

    all_ids = set()
    for ids in results.values():
        all_ids |= ids

    replacements = generate_replacements(all_ids)

    for id in sorted(all_ids):
        print '%s: %s' % (id, replacements[id])

    #for file in sorted(results.keys()):
    #    for id in results[file]:
    #        print '%s: %s' % (file, id)

    return 0


if __name__ == '__main__':
    sys.exit(main())
