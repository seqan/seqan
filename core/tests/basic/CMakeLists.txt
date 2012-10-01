cmake_minimum_required (VERSION 2.6)
project (core_tests_basic)

seqan_add_test_executable(
  test_basic_metaprogramming
  test_basic_metaprogramming.cpp
  test_basic_metaprogramming_logic.h
  test_basic_metaprogramming_type.h)

seqan_add_test_executable(
  test_basic_fundamental
  test_basic_fundamental.cpp
  test_basic_fundamental_helpers.h
  test_basic_fundamental_tags.h
  test_basic_fundamental_metafunctions.h
  test_basic_fundamental_transport.h
  test_basic_fundamental_comparison.h
  test_basic_fundamental_conversion.h
  test_basic_array_construct_destruct.h
  test_basic_hosted_type_interface.h)

seqan_add_test_executable(
  test_basic_concept
  test_basic_concept.cpp
  test_basic_fundamental_concepts.h)

seqan_add_test_executable(
  test_basic_alphabet
  test_basic_alphabet.cpp
  test_basic_alphabet_adapt_builtins.h
  test_basic_alphabet_bio.h
  test_basic_alphabet_concepts.h
  test_basic_alphabet_math.h
  test_basic_alphabet_profile.h
  test_basic_alphabet_qualities.h
  test_basic_alphabet_residue.h
  test_basic_alphabet_storage.h
  test_basic_fundamental_helpers.h)

seqan_add_test_executable(
  test_basic_aggregate
  test_basic_aggregate.cpp
  test_basic_aggregate.h)

seqan_add_test_executable(
  test_basic_allocator
  test_basic_allocator.cpp
  test_basic_allocator.h)

seqan_add_test_executable(
  test_basic_parallelism
  test_basic_parallelism.cpp
  test_basic_parallelism.h)

seqan_add_test_executable(
  test_basic_math
  test_basic_math.cpp
  test_basic_math.h)

seqan_add_test_executable(
  test_basic_smart_pointer
  test_basic_smart_pointer.cpp
  test_basic_smart_pointer_holder.h)

seqan_add_test_executable(
  test_basic_container
  test_basic_container.cpp)

seqan_add_test_executable(
  test_basic_proxy
  test_basic_proxy.cpp
  test_basic_proxy.h)

seqan_add_test_executable(
  test_basic_iterator
  test_basic_iterator.cpp
  test_basic_iterator.h)
