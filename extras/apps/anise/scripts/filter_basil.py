#!/usr/bin/env python
"""Filtration of BASIL output VCF.

Usage: filter_basil.py [OPTIONS] -i basil.vcf -o basil.filtered.vcf
"""

from __future__ import print_function

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'


import argparse
import sys


class BasilVcfFilterConf(object):
    def __init__(self, min_oea_each, min_oea_one, min_oea_sum,
                 min_clipping_each, min_clipping_sum, min_gscore):
        self.min_oea_each = min_oea_each
        self.min_oea_one = min_oea_one
        self.min_oea_sum = min_oea_sum
        self.min_clipping_each = min_clipping_each
        self.min_clipping_sum = min_clipping_sum
        self.min_gscore = min_gscore

    def printConf(self):
        print('min_oea_each\t%s' % self.min_oea_each)
        print('min_oea_one\t%s' % self.min_oea_one)
        print('min_oea_sum\t%s' % self.min_oea_sum)
        print('min_clipping_each\t%s' % self.min_clipping_each)
        print('min_clipping_sum\t%s' % self.min_clipping_sum)
        print('min_gscore\t%s' % self.min_gscore)


class BasilVcfFilter(object):
    def __init__(self, in_file_name, out_file_name, conf):
        self.in_file_name = in_file_name
        self.out_file_name = out_file_name
        self.conf = conf
        self.current_line = None
        self.lines_out = 0

    def run(self):
        with open(self.in_file_name, 'rb') as fin:
            with open(self.out_file_name, 'wb') as fout:
                self._copyHeader(fin, fout)
                self._filter(fin, fout)
        print('Wrote %d lines to output.' % self.lines_out, file=sys.stderr)

    def _copyHeader(self, fin, fout):
        self.current_line = fin.readline()
        while self.current_line and self.current_line.startswith('#'):
            fout.write(self.current_line)
            if 'fileformat' in self.current_line:
                # Write out info about the filter.
                fout.write('##basil_filter=%s\n' % repr(' '.join(sys.argv)))
            self.current_line = fin.readline()

    def _filter(self, fin, fout):
        while self.current_line:
            arr = map(float, self.current_line.split('\t')[9].strip().split(':'))
            gsum, cl_left, cl_right, oea_left, oea_right = arr
            skip = False
            if self.conf.min_gscore:
                if gscore < self.conf.gscore:
                    skip = True
            if self.conf.min_clipping_each:
                if cl_left < self.conf.min_clipping_each or cl_right < self.conf.min_clipping_each:
                    skip = True
            if self.conf.min_clipping_sum:
                if cl_left + cl_right < self.conf.min_clipping_sum:
                    skip = True
            if self.conf.min_oea_each:
                if oea_left < self.conf.min_oea_each or oea_right < self.conf.min_oea_each:
                    skip = True
            if self.conf.min_oea_one:
                if oea_left < self.conf.min_oea_one and oea_right < self.conf.min_oea_one:
                    skip = True
            if self.conf.min_oea_sum:
                if oea_left + oea_right < self.conf.min_oea_sum:
                    skip = True
            if not skip:
                fout.write(self.current_line)
                self.lines_out += 1
            self.current_line = fin.readline()  # read next


def main(args):
    parser = argparse.ArgumentParser(description='Filter BASIL output VCF.')
    parser.add_argument('-i', dest='in_file_name', required=True, help='Input file name.')
    parser.add_argument('-o', dest='out_file_name', required=True, help='Output file name.')

    parser.add_argument('--min-oea-each-side', dest='min_oea_each', type=int,
                        help='Minimal OEA coverage on each side.')
    parser.add_argument('--min-oea-one-side', dest='min_oea_one', type=int,
                        help='Minimal OEA coverage on one side.')
    parser.add_argument('--min-oea-sum', dest='min_oea_sum', type=int,
                        help='Minimal total OEA coverage.')
    parser.add_argument('--min-clipping-each-side', dest='min_clipping_each', type=int,
                        help='Minimal OEA coverage on each side.')
    parser.add_argument('--min-clipping-sum', dest='min_clipping_sum', type=int,
                        help='Minimal total OEA coverage.')
    parser.add_argument('--min-gscore', dest='min_gscore', type=float,
                        help='Smallest geometric mean score')
    args = parser.parse_args()

    conf = BasilVcfFilterConf(min_oea_each=args.min_oea_each,
                              min_oea_sum=args.min_oea_sum,
                              min_oea_one=args.min_oea_one,
                              min_clipping_each=args.min_clipping_each,
                              min_clipping_sum=args.min_clipping_sum,
                              min_gscore=args.min_gscore)
    conf.printConf()

    the_filter = BasilVcfFilter(args.in_file_name, args.out_file_name, conf)
    the_filter.run()
    return 0


if __name__ == '__main__':
  sys.exit(main(sys.argv))

