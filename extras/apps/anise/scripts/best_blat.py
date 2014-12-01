#!/usr/bin/env python
"""Filter BLAT file and obtain best BLAT match for each reference/query.

USAGE: best_blat.py -b in.psl --best-for query_name >out.tsv
"""

from __future__ import print_function

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import argparse
import math
import sys


class PslBlock(object):
    """A block in a BLAT match."""

    def __init__(self, block_size, block_start, target_start):
        self.block_size = block_size
        self.block_start = block_start
        self.target_start = target_start


class PslRecord(object):
    """A PSL entry."""

    def __init__(self, matches, mismatches, rep_matches, ns, query_gap_count,
                 query_gap_bases, target_gap_count, target_gap_bases, strand,
                 query_name, query_size, query_start, query_end,
                 target_name, target_size, target_start, target_end,
                 block_count, block_sizes, query_starts,
                 target_starts):
        self.matches = int(matches)
        self.mismatches = int(mismatches)
        self.rep_matches = int(rep_matches)
        self.ns = int(ns)
        self.query_gap_count = int(query_gap_count)
        self.query_gap_bases = int(query_gap_bases)
        self.target_gap_count = int(target_gap_count)
        self.target_gap_bases = int(target_gap_bases)
        self.strand = strand
        self.query_name = query_name
        self.query_size = int(query_size)
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.target_name = target_name
        self.target_size = int(target_size)
        self.target_start = int(target_start)
        self.target_end = int(target_end)
        self.block_count = int(block_count)
        def split(s):
            return map(int, s.split(',')[:-1])
        self.block_sizes = split(block_sizes)
        self.query_starts = split(query_starts)
        self.target_starts = split(target_starts)
        self.blocks = self.build_blocks(self.block_sizes, self.query_starts, self.target_starts)
        self.compute_stats()

    @classmethod
    def header(klass):
        """Header names."""
        return ['matches', 'mismatches', 'rep_matches', 'ns', 'query_gap_count',
                'query_gap_bases', 'target_gap_count', 'target_gap_bases',
                'strand', 'query_name', 'query_size', 'query_start',
                'query_end', 'target_name', 'target_size', 'target_start',
                'target_end', 'block_count', 'block_sizes', 'query_starts',
                'target_starts']

    def toPsl(self):
        values = [self.matches, self.mismatches, self.rep_matches, self.ns,
                  self.query_gap_count, self.query_gap_bases,
                  self.target_gap_count, self.target_gap_bases,
                  self.strand,
                  self.query_name, self.query_size, self.query_start,
                  self.query_end, self.target_name, self.target_size,
                  self.target_start, self.target_end, self.block_count,
                  ','.join(map(str, self.block_sizes)) + ',',
                  ','.join(map(str, self.query_starts)) + ',',
                  ','.join(map(str, self.target_starts)) + ',']
        return '\t'.join(map(str, values))

    def __repr__(self):
        keys = ['matches', 'mismatches', 'rep_matches',
                'query_name', 'query_start', 'query_end',
                'target_name', 'target_start', 'target_end',
                'query_gap_count', 'target_gap_count', 'block_uncovered']
        tpl = 'PslRecord(' + ', '.join(['%s=%%(%s)s' % (key, key) for key in keys]) + ')'
        vals = dict((key, getattr(self, key)) for key in keys)
        return tpl % vals

    def build_blocks(self, block_sizes, query_starts, target_starts):
        return [PslBlock(*arr) for arr in zip(block_sizes, query_starts, target_starts)]

    def compute_stats(self):
        """Compute number of gaps smaller than 5, count as errors. """
        # TODO(holtgrew): Count leading and trailing query gaps as errors, should not occur.
        self.block_uncovered = 0  # number of uncovered bases in target
        self.block_gap_bases = 0  # gaps in query < 5 and gaps in target of any size
        self.block_mismatch_bases = 0  # uncovered in both, count as mismatch
        self.block_uncovered += self.blocks[0].block_start
        self.block_uncovered += (self.query_size - self.blocks[-1].block_start -
                                 self.blocks[-1].block_size)
        prev_block = self.blocks[0]
        for block in self.blocks[1:]:
            prev_block_query_end = prev_block.block_start + prev_block.block_size
            prev_block_target_end = prev_block.target_start + prev_block.block_size
            target_gaps = block.block_start - prev_block_query_end  # bases in target aligning against gap
            query_gaps = block.target_start - prev_block_target_end  # bases in query aligning against gap
            pseudo_gaps = min(query_gaps, target_gaps)  # not covered by match in both
            #print('target_name', self.target_name,
            #      'query_name', self.query_name,
            #      'target_start', self.target_start,
            #      'target_end', self.target_end,
            #      'pseudo_gaps', pseudo_gaps)
            # Update stats.
            self.block_uncovered += target_gaps
            if pseudo_gaps < 5:  # somehow both not aligned
                self.block_mismatch_bases += pseudo_gaps
                self.block_gap_bases += query_gaps - pseudo_gaps
            # Use this block as previous in next iteration.
            prev_block = block
        # The number of bases in blocks, matches and mismatches.
        self.block_aligned_bases = sum(block.block_size for block in self.blocks) + self.block_gap_bases

    @property
    def length(self):
        """Length without big gaps."""
        return self.matches + self.mismatches + self.rep_matches
        #return self.matches + self.mismatches

    @property
    def score(self):
        return ((self.matches + self.rep_matches / 2) -
                self.mismatches - self.query_gap_count - self.target_gap_count)

    @property
    def coverage(self):
        """Query coverage."""
        return 1.0 * (self.matches + self.mismatches + self.rep_matches) / self.query_size

    @property
    def target_coverage(self):
        """Target coverage."""
        return 1.0 * (self.matches + self.mismatches + self.rep_matches) / self.target_size

    @property
    def identity(self):
        """Compute identity, same as BLAT."""
        q_ali_size = self.query_end - self.query_start
        t_ali_size = self.target_end - self.target_start
        ali_size = min(q_ali_size, t_ali_size)
        if ali_size <= 0:
            return 0.0
        size_dif = q_ali_size - t_ali_size
        if size_dif < 0:
            size_dif = 0
        insert_factor = self.query_gap_count
        total = self.matches + self.rep_matches + self.mismatches
        milli_bad = (1000.0 * (self.mismatches + insert_factor + round(3 * math.log(1 + size_dif)))) / total
        return (100.0 - milli_bad * 0.1) / 100.0
        #return 1.0 * self.matches / (self.matches + self.mismatches)

    @property
    def coords(self):
        """Returns the anchoring genomic location, ready for samtools faidx."""
        return '%s:%d-%d' % (self.target_name, self.target_start + 1, self.target_end)


def loadPsl(fp):
    """Load PSL, yielding PslRecord objects."""
    # Skip header.
    for i in range(0, 5):
        fp.readline()
    for line in fp:
        record = PslRecord(*line.split())
        if record.target_gap_bases > 50000:
            record.target_gap_bases = 0
        yield record


class AniseRecord(object):
    def __init__(self, meta, seq_id, ref, pos, steps, anchored_left,
                 anchored_right, spanning, seq):
        self.meta = meta
        self.seq_id = seq_id
        self.ref = ref
        self.pos = int(pos)
        self.steps = int(steps)
        yesno = {'yes' : True, 'no': False }
        self.anchored_left = yesno[anchored_left]
        self.anchored_right = yesno[anchored_right]
        self.spanning = yesno[spanning]
        self.seq = seq
        self.orig_anchor = None
        self.anchor = None  # anchoring PslRecord


def run(args):
    # Factorize PSL records by query name.
    records = list(loadPsl(args.blat_in_file))
    if args.best_for:
        best = {}
        for record in records:
            if record.identity < args.min_identity:
                continue
            if (not best.get(getattr(record, args.best_for)) or
                best[getattr(record, args.best_for)].score < record.score):
              best[getattr(record, args.best_for)] = record
        records = best.values()
    keys = ['identity', 'query_coverage', 'target_coverage',
            'blat_score', 'novel_bases'] + PslRecord.header()
    print('#' + '\t'.join(keys))
    for record in records:
        print('%.1f\t%.1f\t%.1f\t%d\t%d\t%s' % (100.0 * record.identity,
                                                100.0 * record.coverage,
                                                100.0 * record.target_coverage,
                                                record.score,
                                                record.novel_bases,
                                                record.toPsl()))


def main():
    parser = argparse.ArgumentParser(description='Compute statistics from BLASTN results.')
    parser.add_argument('-b', '--blat-in-file', help='BLAT input file.',
                        type=argparse.FileType('rb'), required=True)
    parser.add_argument('-i', '--min-identity', help='Smallest BLAT identity.',
                        type=float, default=.95)
    parser.add_argument('--best-for', type=str, choices=['query_name', 'target_name'])
    args = parser.parse_args()

    run(args)

    return 0


if __name__ == '__main__':
    sys.exit(main())
