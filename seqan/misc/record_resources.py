#!/usr/bin/env python

import os.path
import sys
# Make Python stuff in py_lib available.
py_lib_path = os.path.join(os.path.dirname(__file__), 'py_lib')
sys.path.append(py_lib_path)
py_lib_path = os.path.join(os.path.dirname(__file__), '..', 'projects', 'benchmarks', 'read_mappers', 'py_lib')
sys.path.append(py_lib_path)

import rmbench.commands as c

USAGE = 'USAGE: record_usage TARGET_FILE PROGRAM_COMMAND_LINE'

def main(args):
    if len(args) <= 2:
        print >>sys.stderr, 'Wrong number of parameters!'
        print >>sys.stderr, USAGE
        return 1

    target_file = args[1]
    binary = args[2]
    arguments = args[3:]

    print 'Executing', binary, ' '.join(arguments)
    print 'target file:', target_file

    c.Task(commands=[c.Command(binary, arguments)], result=target_file).execute()

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
