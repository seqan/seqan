#!/usr/bin/env awk -f

# awk script to get breakpoint contigs
#
# The script works on a breakpoint TSV file as written by mason_variator
# to the path given by the --out-breakpoints.  It then generates a command
# that uses "samtools faidx" to get the infixes around the positions from
# the variation FASTA file.
#
# USAGE: awk -f breakpoint_contigs.awk -v ref=GENOME.fa IN.tsv
#
# Optionally, the variable context can be set to increase the length
# of the context to cut out around breakpoints.

# Set default value for context, initialize cmd in which we will build
# the command.
BEGIN {
    if (context == "")
        context = 100
    cmd = "samtools faidx " ref " ";
}

# For each non-comment record, extend the command for samtools.
/^[^#]/ {
    from = ($3 > context) ? ($3 - context) : 0
    cmd = cmd sprintf(" \"%s:%d-%d\" ", $1, from, $3 + context)
}

# After parsing the file, execute samtools faidx.
END {
    print "EXECUTING " cmd >"/dev/stderr"
    system(cmd)
}
