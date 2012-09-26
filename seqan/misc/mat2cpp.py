#!/usr/bin/env python2.5
"""Simple converter for matrix files to C++ fragments.

A matrix file is read from stdin and appropriate C++ code for
score_matrix_data.h is written to stdout.
"""

import sys


# String for one level of indentation.n
INDENT = "    ";
# Start level of indentation.
INDENT_LEVEL = 2;


def main():
  comments = []
  data_lines = []
  first_data = True
  for line in sys.stdin:
    if line[0] == '#':
      comments.append('// ' + line[1:].strip())
      continue
    if first_data:
      first_data = False
      continue
    data = line.strip().split()[1:]
    formatted_data = ['%3d' % int(d) for d in data]
    data_lines.append(INDENT + ', '.join(formatted_data) + ',')
    

  print '\n'.join([INDENT_LEVEL * INDENT + l for l in comments])
  print INDENT_LEVEL * INDENT + 'static int const _data[TAB_SIZE] = {'
  print '\n'.join([INDENT_LEVEL * INDENT + l for l in data_lines])
  print INDENT_LEVEL * INDENT + '};'
  return 0


if __name__ == '__main__':
  sys.exit(main())
