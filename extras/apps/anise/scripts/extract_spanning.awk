#!/usr/bin/awk -f

# Filter ANISE output FASTA files.
#
# USAGE: filter_spanning.awk <anise.fa >anise.filtered.fa

BEGIN {
    ignore = 0;
    name = "";
    seq = "";
    if (!N_FILTER)
        N_FILTER = "NNNNNN";
    if (!MIN_LEN)
        MIN_LEN = 400;
}

/^>.*SPANNING=yes/ {
    if (!ignore && seq && length(seq) > MIN_LEN && !match(seq, N_FILTER))
        printf(">%s\n%s\n", name, seq);
    ignore = !!match($0, "STOPPED=too_many_reads");
    name = substr($0, 2);
    seq = "";
}

/^>.*SPANNING=no/ {
    if (!ignore && seq && length(seq) > MIN_LEN && !match(seq, N_FILTER))
        printf(">%s\n%s\n", name, seq);
    ignore = 1;
    name = substr($0, 2);
    seq = "";
}

!/^>/ {
    seq = seq $0;
}

END {
    if (!ignore && seq && length(seq) > MIN_LEN && !match(seq, N_FILTER))
        printf(">%s\n%s\n", name, seq);
}
