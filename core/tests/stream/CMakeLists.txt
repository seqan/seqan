cmake_minimum_required (VERSION 2.6)
project (extras_tests_stream)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_stream
                          test_stream.cpp
                          test_stream_adapt_cstdio.h
                          test_stream_adapt_fstream.h
                          test_stream_adapt_sstream.h
                          test_stream_bz2_file.h
                          test_stream_char_array.h
                          test_stream_gz_file.h
                          test_stream_tokenize.h
                          test_stream_generic.h
                          test_stream_lexical_cast.h
                          test_stream_bgzf.h)
