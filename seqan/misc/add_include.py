#!/usr/bin/env python2.5
"""DDDoc Helper Script -- add includes.

Called the path to the SeqAn checkout, this tool will search for all
SeqAn header files.  Within each file, it will look for all
Metafunction, Function, Tag, Class and Specialization dddoc
documentation.  For each of documentation comment, it will add the
line "..include:<seqan/${module}.h>" where ${module} is the name of the module
that defines the symbol.

Note that the ddoc search is purely heuristic.  It will look for the
first line containing "*/" to find the end of a comment.
"""

from __future__ import with_statement

import os
import os.path
import sys

USAGE_STR = """
DDDoc Helper Script -- add includes.

USAGE: add_include.py PATH_TO_SEQAN.
""".strip()


def processHeaderFile(base_path, relative_path):
  """Process the header and insert the .include line if it is not already there."""
  full_path = os.path.join(base_path, relative_path)
  print full_path
  # Read file.
  with open(full_path, 'r') as f:
    contents = f.readlines()

  # Try to detect the line ending.
  if contents[0][-2:] == '\r\n':
    line_ending = '\r\n'
  else:
    assert contents[0][-1] == '\n'
    line_ending = '\n'

  # Now, try loop over the lines and look for a line beginning with
  # .Function, .Metafunction, .Tag, .Class or .Spec.  Setting includes
  # for member variables or functions does not make sense.
  hit = False  # Flag: Changed file contents?
  result = []  # Output file lines.
  in_doc = False  # State: In documentation comment?
  has_include = False  # State: Have we seen a .include line?
  for line in contents:
    if in_doc:
      if line.startswith('..include:'):
        has_include = True
      elif '*/' in line:
        # End of the comment block, maybe insert '.include' line.
        if not has_include:
          hit = True
          module_name = relative_path.split('/')[0]
          result.append('..include:seqan/%s.h%s' % (module_name, line_ending))
        # Reset state automaton.
        in_doc = False
        has_include = False
    else:
      if line.startswith('.Function') or line.startswith('.Metafunction') or \
            line.startswith('.Tag') or line.startswith('.Spec') or \
            line.startswith('.Class'):
        in_doc = True
    result.append(line)

  # Write file again:
  if hit:
    with open(full_path, 'w') as f:
      f.write(''.join(result))


def collectHeaderFilenames(base_path):
  """Return all header files below the given base path in a list."""
  result = []
  for root, dirs, files in os.walk(base_path):
    for f in files:
      if not f.endswith('.h'):
        continue  # Ignore non-headers.
      if f.startswith('.'):
        continue  # Ignore hidden files
      result.append(os.path.join(root, f)[len(base_path) + 1:])
    # Do not visit .svn directories.
    if '.svn' in dirs:
      dirs.remove('.svn')
  return result


def main():
  if len(sys.argv) != 2:
    print >>sys.stderr, "ERROR: Wrong argument count."""
    print >>sys.stderr, USAGE_STR
    return 1

  base_path = sys.argv[1]
  if base_path[-1] == '/':
    base_path = base_path[-1]
  header_rel_paths = collectHeaderFilenames(base_path)

  for header_rel_path in header_rel_paths:
    print header_rel_path
    processHeaderFile(base_path, header_rel_path)

if __name__ == '__main__':
  sys.exit(main())
