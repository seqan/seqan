#!/usr/bin/env python2
"""Generation of detailed ROI reports with larger plots.

This report generation works for hundred of ROIs.
"""

try:
    import argparse
except ImportError:
    import argparse26 as argparse
import math
import os.path
import sys

import Cheetah.Template
import matplotlib.pyplot as plt

import ngs_roi.app
import ngs_roi.argparse
import ngs_roi.io

PAGE_TPL = """
<html>
  <head>
    <title>ROI Table</title>
    <style type="text/css">
    div.plot
    {
        float: left;
        padding: 4px;
        margin: 2px;
        width: 420px;
    }

    .plot h2 { margin-top: 3px; margin-bottom: 3px; text-align: center; }
    .plot img { display: block; margin: 0 auto; }
    </style>
  </head>
  <body>
    <h1>Detailed ROI Report</h1>

    #for i, roi in enumerate($records)
    <div class="plot">
      <h2>${roi.ref}:${roi.start_pos + 1}-${roi.end_pos+1}</h2>
      <a href="${href($roi)}" target="dead"><img src="plot_${i}.png" /></a>
      <p>
        <b>chr:start-end</b> <a href="${href($roi)}" target="dead">${roi.ref}:${roi.start_pos}-${roi.end_pos} ${roi.strand}</a>;
        <b>region name</b> ${roi.region_name};
        <b>region length</b> ${roi.region_length};
      </p>
      #if $roi.data
      <p>#for j, key in enumerate($data_keys)#<b>$key:</b> ${roi.data[$j]}; #end for#</p>
      #end if
    </div>
    #end for
    <iframe name="dead" height="0" width="0"></iframe>
    <div><code>$args</code></div>
  </body>
</html>
"""

class DetailedRoiGenerator(ngs_roi.app.App):
    """Generate detailed ROI report.

    :ivar args:Arguments from the comment line.
    """

    def __init__(self, args):
        self.args = args

    def run(self):
        """Run report generation, return status code.

        :return: integer with the result.
        """
        print >>sys.stderr, 'Loading ROI'
        records = ngs_roi.io.load(self.args.in_file, self.args.max_rois)
        keys = records[0].data_keys

        self.writeHtml(keys, records)
        self.writePlots(records)
        return 0

    def writePlots(self, records):
        COLOR = 'blue'
        LINE_WIDTH = .5
        LINE_STYLE = '-'
        TICK_FONT_SIZE = 8
        LABEL_FONT_SIZE = 10
        for i, roi in enumerate(records):
            file_name = 'plot_%d.png' % i
            file_name = os.path.join(self.args.out_dir, file_name)
            print >>sys.stderr, 'Writing plot %s' % file_name
            plt.figure(figsize=(4, 2.5))
            plt.gcf().subplots_adjust(bottom=0.16, left=0.15)
            plt.plot(roi.points, color=COLOR, linewidth=LINE_WIDTH, linestyle=LINE_STYLE)
            plt.ylim(ymin=0)
            if self.args.max_value:
                plt.ylim(ymax=self.args.max_value)
            plt.tick_params(labelsize=TICK_FONT_SIZE)
            plt.ylabel('coverage', fontsize=LABEL_FONT_SIZE, weight='semibold')
            plt.xlabel('ROI beginPos', fontsize=LABEL_FONT_SIZE, weight='semibold')
            plt.savefig(file_name)

    def writeHtml(self, keys, records):
        file_name = self.args.out_file
        print >>sys.stderr, 'Writing HTML file %s' % file_name

        vals = {'args': self.args, 'records': records, 'data_keys': keys,
                'href': lambda x: self.buildHref(x.ref, x.start_pos, x.end_pos)}
        t = Cheetah.Template.Template(PAGE_TPL, searchList=vals)
        
        with open(file_name, 'wb') as f:
            f.write(str(t))


def main():
    parser = argparse.ArgumentParser(description='Plot ROI file.')
    ngs_roi.argparse.addFileArguments(parser)
    ngs_roi.argparse.addPlotGridArguments(parser)
    ngs_roi.argparse.addLinkArguments(parser)
    args = parser.parse_args()
    ngs_roi.argparse.applyFileDefaults(args)

    app = DetailedRoiGenerator(args)
    return app.run()


if __name__ == '__main__':
    sys.exit(main())
