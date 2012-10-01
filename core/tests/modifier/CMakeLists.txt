cmake_minimum_required (VERSION 2.6)
project (core_tests_modifier)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_modifier test_modifier.cpp
					helpers.h
					test_modifier_alphabet.h
					test_modifier_functors.h
					test_modifier_shortcuts.h
					test_modifier_view.h)
