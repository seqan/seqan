cmake_minimum_required (VERSION 2.6)
project (extras_tests_bam_io)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_bam_io
                          test_bam_io.cpp
                          test_bam_alignment_record.h
                          test_bam_header_record.h
                          test_bam_io_context.h
                          test_bam_sam_conversion.h
                          test_bam_tags_dict.h
                          test_read_sam.h
                          test_write_bam.h
                          test_write_sam.h
                          test_bam_stream.h)
