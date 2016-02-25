#!/usr/bin/env python2
"""Create ROI overview report.

This report consists of plots of all metrics (y: metric, x: rank of value).
Each plot is written out as a PNG file and we also create one output HTML file
that shows all HTML files.

Plotting is done using the fine matplotlib.
"""

from __future__ import print_function

__author__    = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'
__copyright__ = 'Copyring 2013, Freie Universitaet Berlin'
__license__   = 'BSD 3-clause'


import sys
import os
import os.path
try:
    import argparse
except ImportError:
    import argparse26 as argparse

import Cheetah.Template
import matplotlib.pyplot as plt

import ngs_roi.io
import ngs_roi.app
import ngs_roi.argparse


# The HTML template to use for generating the HTML page.
REPORT_TPL = """
<html>
  <head>
    <title>ROI Report</title>
  </head>
  <body>
    <h1>ROI Report</h1>
    <h2>Table of Contents</h2>
    <ul>
    #for figure in $figures
      <li><a href="#$figure.slug">$figure.title</a></li>
    #end for
    </ul>
    <h2>Plots</h2>
    #for figure in $figures
    <h3 id="$figure.slug">$figure.title</h3>
    <img src="$figure.file_name" title="$figure.title" />
    #end for
  </body>
</html>
"""


class ReportBuilder(ngs_roi.app.App):
    """This class is used for building the report."""

    def __init__(self, args):
        self.args = args
        self.in_file = self.args.in_file
        self.out_file = self.args.out_file
        self.out_dir = self.args.out_dir
        self.prepareOutDir()

    def plotAndWrite(self, file_name, numbers, ylabel):
        """Create plot of numbers and write as PNG.

        :param file_name:path to write plot to as image
        :param numbers:list of numbers to plot
        :param ylabel:label for y axis
        """
        plt.figure()
        plt.plot(numbers)
        plt.ylabel(ylabel)
        plt.savefig(file_name)

    def run(self):
        # Load ROI.
        print('Loading ROI...', file=sys.stderr)
        records = ngs_roi.io.load(self.in_file)
        keys = records[0].data_keys
        # Create ROI plots.
        print('Creating plots...', file=sys.stderr)
        METRICS = [('start position', 'start_pos', lambda x: x.start_pos),
                   ('end position', 'end_pos', lambda x: x.end_pos),
                   ('region length', 'region_length', lambda x: x.region_length),
                   ('max count', 'max_count', lambda x: x.max_count)]
        def getData(i):
            def func(x):
                try:
                    res = float(x.data[i])
                except ValueError:
                    res = 0
                return res
            return func
        for i, key in enumerate(keys):
            slug = ''.join(x for x in key if x.isalnum())
            METRICS.append((key, slug, getData(i)))
        figure_infos = []
        for title, slug, func in METRICS:
            values = [func(x) for x in records]
            file_name = 'report_%s.png' % slug
            file_path = os.path.join(self.out_dir, file_name)
            self.plotAndWrite(file_path, sorted(values), title)
            figure_infos.append({'file_name': file_name, 'title': title, 'slug': slug})
        # Create report HTML.
        name_space = {'figures': figure_infos}
        t = Cheetah.Template.Template(REPORT_TPL, searchList=name_space)
        with open(os.path.join(self.out_dir, 'index.html'), 'wb') as f:
            f.write(str(t))
        with open(os.path.join(self.out_file), 'wb') as f:
            f.write(str(t))

def main():
    """Program entry point."""

    parser = argparse.ArgumentParser(description='Create ROI report.')

    ngs_roi.argparse.addFileArguments(parser)
    args = parser.parse_args()
    ngs_roi.argparse.applyFileDefaults(args)

    report_builder = ReportBuilder(args)
    return report_builder.run()


if __name__ == '__main__':
    sys.exit(main())