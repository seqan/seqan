
# Hacky... SAMTOOLS_FOUND set in top directory.
if (ZLIB_FOUND)
  SET(LIBBAM_SRC bgzf.c
                 kstring.c
                 bam_aux.c
                 bam.c
                 bam_import.c
                 sam.c
                 bam_index.c
                 bam_pileup.c
                 bam_lpileup.c
                 bam_md.c
                 glf.c
                 razf.c
                 faidx.c
                 knetfile.c
                 bam_sort.c
                 sam_header.c
                 bam_reheader.c
                 kprobaln.c)

  ADD_DEFINITIONS(-D_USE_KNETFILE -D_FILE_OFFSET_BITS=64)
  INCLUDE_DIRECTORIES(${SEQAN_LIBRARY_ABS}/lib/samtools)
  ADD_LIBRARY(bam ${LIBBAM_SRC})
else (ZLIB_FOUND)
  message("WARNING: No ZLIB: samtools disabled")
endif (ZLIB_FOUND)
