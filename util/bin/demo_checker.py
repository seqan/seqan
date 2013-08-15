#!/usr/bin/env python
"""Demo checker script.

Given a demo .cpp file PATH.cpp we can make it a small test if there is a file
PATH.cpp.stdout and/or PATH.cpp.stderr.  The test is implemented using this
script.

The script is called with the options --binary-path and one or both of
--stdout-path and --stderr-path.  The demo is executed and the test succeeds
if the exit code is 0 and the standard/error output is the same as in the
.stdout/.stderr file.  If there is output and the file is missing then this is
a failure as well.
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'


import argparse
import difflib
import subprocess
import sys


def t(s):
    """Force Windows line endings to Unix line endings."""
    return s.replace("\r\n", "\n")


def loadExpected(args):
    """Load the expected file contents."""
    out, err = '', ''
    if args.stdout_path:
        with open(args.stdout_path, 'r') as f:
            out = f.read()
    if args.stderr_path:
        with open(args.stderr_path, 'r') as f:
            err = f.read()
    return t(out), t(err)


def runDemo(args):
    cmd = [args.binary_path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    return t(p.stdout.read()), t(p.stderr.read()), p.returncode


def main():
    """Program entry point."""
    parser = argparse.ArgumentParser(description='Run SeqAn demos as apps.')
    parser.add_argument('--binary-path', dest='binary_path', required='True',
                        help='Path to the demo binary to execute.')
    parser.add_argument('--stdout-path', dest='stdout_path',
                        help='Path to standard out file to compare to.',
                        default=None)
    parser.add_argument('--stderr-path', dest='stderr_path',
                        help='Path to standard error file to compare to.',
                        default=None)
    args = parser.parse_args()

    print >>sys.stderr, 'Loading files "%s", "%s".' % (args.stdout_path, args.stderr_path)
    expected_out, expected_err = loadExpected(args)
    print >>sys.stderr, 'Running %s.' % args.binary_path
    actual_out, actual_err, ret = runDemo(args)

    if ret != 0:
        print >>sys.stderr, 'ERROR: Return code of %s was %s.' % (args.binary_path, ret)
        return 1
    else:
        print >>sys.stderr, 'Return code was %s.' % ret

    if expected_out != actual_out:
        print >>sys.stderr, 'The standard output was not as expected!'
        l = difflib.context_diff(expected_out, actual_out,
                                 fromfile='expected', tofile='actual')
        print >>sys.stderr, '\n'.join(l)
    else:
        print >>sys.stderr, 'Standard output was as expected.'

    if expected_err != actual_err:
        print >>sys.stderr, 'The standard errput was not as expected!'
        l = difflib.context_diff(expected_err, actual_err,
                                 fromfile='expected', tofile='actual')
        print >>sys.stderr, '\n'.join(l)
    else:
        print >>sys.stderr, 'Standard error was as expected.'

    return not (expected_out == expected_out and expected_err == actual_err)


if __name__ == '__main__':
    sys.exit(main())
