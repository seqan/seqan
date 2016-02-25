#!/usr/bin/env python2
"""ROI Table Generator

Generates a HTML page with a table of ROI record details.  Besides showing the
numeric ROI information, it also gives little roi plots in a column.

For the little ROI plots, it calls the program roi_plot_thumbnails that has to
be in the PATH.
"""

__author__    = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'
__copyright__ = 'Copyring 2013, Freie Universitaet Berlin'
__license__   = 'BSD 3-clause'


# TODO(holtgrew): Actually call roi_plot_thumbnails
# TODO(holtgrew): from __future__ use print_function


try:
    import argparse
except ImportError:
    import argparse26 as argparse
import math
import os.path
import sys

import Cheetah.Template

import ngs_roi.app
import ngs_roi.argparse
import ngs_roi.io


# Main template.
PAGE_TPL = """
<html>
  <head><title>ROI Table</title></head>
  <body>
    <h1>ROI Table</h1>
    $table
    <div><code>$args</code></div>
  </body>
</html>
"""

# Template for a table.
TABLE_TPL = """
<table border="1">
  <tr>
    <th>plot</th>
    <th>chr</th>
    <th>start</th>
    <th>end</th>
    <th>name</th>
    <th>length</th>
    <th>strand</th>
    <th>max_count</th>
    #for i, key in enumerate($data_keys)
    <th>$key</th>
    #end for
  </tr>
  #for id, roi in enumerate($records)
  <tr>
    <td><div style="width:${args.plot_width}; margin:2px; height:${args.plot_height+1}; background:url(thumbnail_${imgId($id)}.png) -${imgX($id)} -${imgY($id)};"></div></td>
    <td>$roi.ref</td>
    <td style="text-align:right;">$fmtPos($roi.start_pos + 1)</td>
    <td style="text-align:right;">$fmtPos($roi.end_pos)</td>
    <td><a href="$href($roi)">$roi.region_name</a></td>
    <td style="text-align:right;">$fmtPos($roi.region_length)</td>
    <td style="text-align:center;">$roi.strand</td>
    <td style="text-align:right;">$roi.max_count</td>
    #for i, key in enumerate($data_keys)
    <td>$roi.data[$i]</td>
    #end for
  </tr>
  #end for
</table>
"""


class RoiTable(object):
    """A table of ROI records with small plots."""

    def __init__(self, args, keys, records, app):
        self.args = args
        self.keys = keys
        self.records = records
        self.app = app

    def tplFuncs(self):
        def intWithCommas(x):
            if type(x) not in [type(0), type(0L)]:
                raise TypeError("Parameter must be an integer.")
            if x < 0:
                return '-' + intWithCommas(-x)
            result = ''
            while x >= 1000:
                x, r = divmod(x, 1000)
                result = ",%03d%s" % (r, result)
            return "%d%s" % (x, result)

        def imgId(idx):
            """Image id from roi record id."""
            return idx / (self.args.num_rows * self.args.num_cols)

        def imgX(idx):
            """x position in image from record id."""
            x = idx % self.args.num_cols
            res = x * self.args.plot_width
            if x > 0:
                res += (x - 1) * 2
            return res

        def imgY(idx):
            """y position in image from record id."""
            y = idx / self.args.num_cols
            res = y * self.args.plot_height
            res += y * 2
            return res

        return {'fmtPos': intWithCommas, 'imgId': imgId, 'imgX': imgX, 'imgY': imgY}

    def render(self):
        """Returns string with rendered table."""
        vals = {'data_keys': self.keys, 'records': self.records, 'args': self.args,
                'href': lambda x:self.app.buildHref(x.ref, x.start_pos, x.end_pos)}
        vals.update(self.tplFuncs())
        t = Cheetah.Template.Template(TABLE_TPL, searchList=vals)
        return str(t)


class TableApp(ngs_roi.app.App):
    def __init__(self, args):
        # Call parent's constructor and create output directory.
        ngs_roi.app.App.__init__(self, args)
        self.prepareOutDir()

    def run(self):
        # Load ROI records.
        print >>sys.stderr, 'Loading ROI'
        records = ngs_roi.io.load(self.args.in_file, self.args.max_rois)
        keys = []
        if records:
            keys = records[0].data_keys

        # Create plots.
        print >>sys.stderr, 'Creating plots...'
        runner = ngs_roi.app.PlotThumbnailsRunner(self.args)
        runner.run()

        # Create table.
        print >>sys.stderr, 'Creating table...'
        self.createHtml(self.args.out_file, keys, records)
        return 0

    def createHtml(self, file_name, keys, records):
        print >>sys.stderr, 'Writing %s' % self.args.out_file
        vals = {'table': RoiTable(self.args, keys, records, self).render(),
                'args': self.args}
        t = Cheetah.Template.Template(PAGE_TPL, searchList=vals)
        with open(self.args.out_file, 'wb') as f:
            f.write(str(t))


def main():
    parser = argparse.ArgumentParser(description='Plot ROI file.')

    ngs_roi.argparse.addFileArguments(parser)
    ngs_roi.argparse.addPlotGridArguments(parser, default_plot_height=60,
                                          default_plot_width=90)
    ngs_roi.argparse.addLinkArguments(parser)
    args = parser.parse_args()
    ngs_roi.argparse.applyFileDefaults(args)

    app = TableApp(args)
    return app.run()


if __name__ == '__main__':
    sys.exit(main())
