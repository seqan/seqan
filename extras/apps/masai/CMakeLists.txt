cmake_minimum_required (VERSION 2.6)
project (extras_apps_masai)

OPTION(SEQAN_EXTRAS_MASAI_DISABLE_MMAP "Set to ON to disable memory mapping of the reference genome index." OFF)

IF(SEQAN_EXTRAS_MASAI_DISABLE_MMAP)        
    add_definitions(-DSEQAN_EXTRAS_MASAI_DISABLE_MMAP=1)
ENDIF(SEQAN_EXTRAS_MASAI_DISABLE_MMAP)

seqan_add_executable(masai_indexer indexer.cpp
                                   indexer.h
                                   index.h
                                   store.h
                                   options.h)

seqan_add_executable(masai_mapper mapper.cpp
                                  mapper.h
                                  tags.h
                                  index.h
                                  store.h
                                  indexer.h
                                  seeder.h
                                  extender.h
                                  manager.h
                                  matches.h
                                  stream.h
                                  writer.h
                                  options.h
                                  index/index_qgram_stretched.h
                                  index/find_backtracking_stretched.h
                                  store/store_io.h
                                  store/store_io_sam.h)

seqan_add_executable(masai_output_se sorter.cpp
                                     sorter.h
                                     tags.h
                                     store.h
                                     indexer.h
                                     matches.h
                                     stream.h
                                     writer.h
                                     options.h
                                     store/store_io.h
                                     store/store_io_sam.h)

seqan_add_executable(masai_output_pe pairer.cpp
                                     pairer.h
                                     tags.h
                                     store.h
                                     indexer.h
                                     matches.h
                                     stream.h
                                     writer.h
                                     options.h
                                     store/store_io.h
                                     store/store_io_sam.h)
