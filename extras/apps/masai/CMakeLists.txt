cmake_minimum_required (VERSION 2.6)
project (sandbox_esiragusa_apps_masai)

seqan_add_executable(masai_indexer indexer.cpp indexer.h index.h store.h options.h)
seqan_add_executable(masai_mapper mapper.cpp mapper.h tags.h index.h store.h indexer.h seeder.h extender.h manager.h matches.h writer.h options.h index/index_qgram_stretched.h index/find_backtracking_stretched.h store/store_io.h store/store_io_sam.h)
seqan_add_executable(masai_output_se sorter.cpp sorter.h tags.h store.h indexer.h matches.h writer.h options.h store/store_io.h store/store_io_sam.h)
seqan_add_executable(masai_output_pe pairer.cpp pairer.h tags.h store.h indexer.h matches.h writer.h options.h store/store_io.h store/store_io_sam.h)
