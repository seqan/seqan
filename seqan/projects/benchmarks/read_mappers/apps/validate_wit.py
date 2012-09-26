#!/usr/bin/env python2.5
"""Weighted Interval Target file validator.

Takes an arbitrarily long list of file names and validates them as WIT
files.  The standard input is given as '-'.

Usage:
  validate_wit.py [FILE.wit]+
"""

from __future__ import with_statement

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import sys
import re


def PrintError(filename, line, error):
  """Prints an error message for an error.

  Arg:
    filename  String, path to file.
    line      Line with the errror.
    error     Error message.
  """
  msg = '%(filename)s:%(line)s %(error)s' % {'filename': filename,
                                             'line': line,
                                             'error': error}
  print >>sys.stderr, msg


def ValidateFile(file):
  """Reads in a file and validates it as a WIT file.
  
  Args:
    The file to validate.

  Returns:
    True if the file is valid, False otherwise.
  """
  # Validate header.
  line = file.readline().rstrip()
  header_re = "^@HD\tVN:1.0$"
  if not re.match(header_re, line):
    PrintError(filename, 1, 'Invalid header: "%s".' % line)
    return False
  # Validate the rest of the file.
  lastline = False  # Flag, last line can be empty.
  lineno = 1
  line_re = '^[^ \t\n\r]+\t-?[0-9]+\t[^\t\n\r@=]+\t[0-9]+\t[0-9]+$'
  for line in file:
    lineno += 1
    line = line.rstrip()
    # Line should not be empty if it is not the last line.  Detection
    # is done with the state flag lastline.
    if lastline:
      PrintError(filename, lineno, 'Only the last line can be empty.')
      return False
    if not line:
      lastline = True
    # Skip comments.
    if line[0] == '#':
      continue
    # Validate line.
    if not re.match(line_re, line):
      PrintError(filename, lineno, 'Invalid line: "%s".' % line)
      return False
  return True


def Main(argv):
  """Main entry point.
    
  Args:
    argv  Array with program's arguments.

  Returns:
    0 if everything went fine, an error code != 0 otherwise.
  """
  all_valid = True
  for filename in argv[1:]:
    # Load lines from stdin or the given filename.
    if filename == '-':
      valid = ValidateFile(sys.stdin)
    else:
      with open(filename) as file:
        valid = ValidateFile(file)
    # Validation output for each file.
    if valid:
      print '%s OK' % filename
    all_valid = all_valid and valid
  # Validation output for all files is the return code.
  if all_valid:
    return 0
  else:
    return 1


if __name__ == '__main__':
    sys.exit(Main(sys.argv))
