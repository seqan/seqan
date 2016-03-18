#!/usr/bin/env python2
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

__author__ = """Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
                Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
            """


import argparse
import difflib
import subprocess
import sys
import re

def t(s):
    """Force Windows line endings to Unix line endings."""
    return s.replace("\r\n", "\n")

def fuzzyEqual(pattern, text):
    """checks if the expected output is eqal to the actualoutput using a reqex
        use the literal [VAR] if the part of the output is not expected to be the same all the time.
    """
    if len(pattern) != len(text):
        print >> sys.stderr, 'Number of lines differ. Expected output has %s lines whereas actual has %s lines.' % (len(pattern), len(text))
        return False
    for i in range(len(pattern)):
        T = text[i]
        P = pattern[i]
        if T == P :
            continue
        else :
            if '[VAR]' not in P:
                print >> sys.stderr, 'Line %s is different between expected and actual outputs.' % (i)
                return False
            else:
                P = (re.escape(P)).replace('\\[VAR\\]', "[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?")
                r = re.compile(P)
                if re.match(r, T) == None:
                    print >> sys.stderr, 'Line %s is different (REGEX) between expected and actual outputs.' % (i)
                    return False
    return True

def loadExpected(args):
    """Load the expected file contents."""
    out, err = '', ''
    if args.stdout_path:
        with open(args.stdout_path, 'rb') as f:
            out = f.read()
    if args.stderr_path:
        with open(args.stderr_path, 'rb') as f:
            err = f.read()
    return t(out.strip()).split('\n'), t(err.strip()).split('\n')


def runDemo(args):
    cmd = [args.binary_path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdoutbuff, stderrbuff = p.communicate()
    return t(stdoutbuff.strip()).split('\n'), t(stderrbuff.strip()).split('\n'), p.returncode


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

    print >>sys.stderr, 'Running %s.' % args.binary_path
    actual_out, actual_err, ret = runDemo(args)

    if ret != 0:
        print >>sys.stderr, 'ERROR: Return code of %s was %s.' % (args.binary_path, ret)
        return 1
    else:
        print >>sys.stderr, 'Return code was %s.' % ret

    print >>sys.stderr, 'Loading files "%s", "%s".' % (args.stdout_path, args.stderr_path)
    expected_out, expected_err = loadExpected(args)
    is_stdout_as_expected = fuzzyEqual(expected_out, actual_out)
    is_stderr_as_expected = fuzzyEqual(expected_err, actual_err)

    if not is_stdout_as_expected:
        print >>sys.stderr, 'The standard output was not as expected!'
        l = difflib.context_diff(expected_out, actual_out,
                                 fromfile='expected', tofile='actual')
        print >>sys.stderr, '\n'.join(l)
    else:
        print >>sys.stderr, 'Standard output was as expected.'

    if not is_stderr_as_expected:
        print >>sys.stderr, 'The standard error was not as expected!'
        l = difflib.context_diff(expected_err, actual_err,
                                 fromfile='expected', tofile='actual')
        print >>sys.stderr, '\n'.join(l)
    else:
        print >>sys.stderr, 'Standard error was as expected.'

    # here we used not because we need return-code 0 (False) if test is successful 
    return not (is_stdout_as_expected and is_stderr_as_expected)


if __name__ == '__main__':
    sys.exit(main())
