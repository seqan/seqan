cmake_minimum_required (VERSION 2.6)
project (sandbox_arg_parse_tests_arg_parse)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_arg_parse
                          test_extensions.h
                          test_arg_parse.cpp
                          test_arg_parse.h
                          test_arg_parse_argument.h
                          test_arg_parse_option.h
                          test_argument_parser.h)
