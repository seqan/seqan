#!/usr/bin/env python2
"""Thumbnail Plot Generator.

This report generator uses the binary roi_plot_thumbnails (C++ program, must
be in PATH) to generate many PNG images with small ROI plots.  It then creates
a HTML file that includes the PNG and adds an overlay image link map such that
when the user clicks on a plot, she is transfered to the according position in
a local IGV or the UCSC genome browser, for example.
"""

__author__    = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'
__copyright__ = 'Copyring 2013, Freie Universitaet Berlin'
__license__   = 'BSD 3-clause'

# TODO(holtgrew): Use Cheetah templates to generate HTML.
# TODO(holtgrew): Call roi_plot_thumbnails from Python script.
# TODO(holtgrew): from __future__ use print_function


try:
    import argparse
except ImportError:
    import argparse26 as argparse
import math
import os.path
import sys

import ngs_roi.app
import ngs_roi.argparse
import ngs_roi.io


class LinkRegion(object):
    """Region on picture with genomic interval."""

    def __init__(self, x1, y1, x2, y2, ref, begin_pos, end_pos, name):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.ref = ref
        self.begin_pos = begin_pos
        self.end_pos = end_pos
        self.name = name


class RoiPlotGrid(object):
    """A grid of ROI plots.

    width -- width one one plot
    height -- height of one plot
    """

    # Border in each direction.
    BORDER = 0
    # Spacing between plots.
    SPACE = 2

    def __init__(self, width, height, columns, rows, bg_plot=(0, 0, 0, 10),
                 frame_col=(80, 80, 80, 255), fg_col=(0, 0, 0, 255)):
        # Properties of the grid.
        self.width = float(width)
        self.height = float(height)
        self.columns = columns
        self.rows = rows
        # State.
        self.idx = 0  # currently drawn plot
        self.canvas_width = 2 * self.BORDER + columns * width + (columns - 1) * self.SPACE
        self.canvas_height = 2 * self.BORDER + rows * height + (rows - 1) * self.SPACE
        # Colors.
        self.bg_plot = bg_plot
        self.frame_col = frame_col
        self.fg_color = fg_col
        # List of link region for image map.
        self.link_regions = []

    def plotStart(self, idx):
        """Return pair with start coordinates of plot."""
        row = idx / self.columns
        col = idx % self.columns
        x = self.BORDER + col * (self.width + self.SPACE)
        y = self.BORDER + row * (self.height + self.SPACE)
        return x, y

    def plotRecord(self, record):
        """Register plotting of a record."""
        start_x, start_y = self.plotStart(self.idx)
        self.idx += 1
        # Register link region.
        self.link_regions.append(LinkRegion(
                start_x, start_y, start_x + self.width, start_y + self.height,
                record.ref, record.start_pos, record.end_pos, record.region_name))


class GridLinks(object):
    """Link information for one grid."""

    def __init__(self, file_name, link_regions):
        self.file_name = file_name
        self.link_regions = link_regions


class PlotThumbnailsApp(ngs_roi.app.App):
    def __init__(self, args):
        # Call parent's constructor.
        ngs_roi.app.App.__init__(self, args)
        self.prepareOutDir()
        # Initialize members.
        self.grid = None
        self.plot_idx = 0
        self.grid_links = []

    def run(self):
        # Load ROI records.
        print >>sys.stderr, 'Loading ROI'

        # Create plots.
        runner = ngs_roi.app.PlotThumbnailsRunner(self.args)
        runner.run()

        # Create HTML.
        print >>sys.stderr, 'Creating HTML...'
        num_plots = self.args.num_cols * self.args.num_rows  # plots on grid
        for i, roi in enumerate(ngs_roi.io.RoiFile(self.args.in_file)):
            if self.args.max_rois > 0 and i >= self.args.max_rois:
                break
            # Write out old grid (if any) and create new one.
            if i % num_plots == 0:
                if self.grid:
                    print >>sys.stderr, '  Writing plot %d...' % self.plot_idx
                self.writeGrid()
                self.grid = RoiPlotGrid(self.args.plot_width, self.args.plot_height,
                                        self.args.num_cols, self.args.num_rows)
            # Put the next plot on the grid.
            self.grid.plotRecord(roi)
        print >>sys.stderr, '  Writing plot %d...' % self.plot_idx
        self.writeGrid()  # Write last grid.
        self.createHtml(self.args.out_file)
        return 0

    def writeGrid(self):
        """Register writing of grid."""
        if not self.grid:
            return
        # Append grid info.
        file_name = 'thumbnail_%d.png' % self.plot_idx
        file_name = os.path.join(self.args.out_dir, file_name)
        self.plot_idx += 1
        self.grid_links.append(GridLinks(os.path.basename(file_name), self.grid.link_regions))

    def createHtml(self, file_name):
        print >>sys.stderr, 'Writing HTML to %s' % file_name
        with open(file_name, 'wb') as f:
            f.write('<html><body>\n')
            f.write('<h1>ROI Thumbnail Plots</h1>')
            for gl in self.grid_links:
                vals = (gl.file_name, gl.file_name, self.grid.canvas_width, self.grid.canvas_height)
                f.write('<img src="%s" usemap="#%s" width="%d" height="%d" />\n' % vals)
                f.write('<map name="%s">\n' % gl.file_name)
                for lr in gl.link_regions:
                    locus_label = (lr.ref, lr.begin_pos + 1, lr.end_pos, lr.name)
                    vals = {'x1': lr.x1, 'x2': lr.x2, 'y1': lr.y1, 'y2': lr.y2,
                            'title': '%s %d-%d (%s)' % locus_label,
                            'href': self.buildHref(lr.ref, lr.begin_pos, lr.end_pos),
                            'onclick': ''}
                    # Add onclick handler to prevent opening of new window.
                    vals['target_attr'] = ''
                    if self.args.link_target:
                        vals['target_attr'] = ' target="%s"' % self.args.link_target
                    if self.args.link_type == 'local_igv':
                        vals['target_attr'] = ' target="empty"'
                    f.write('  <area shape="rect" coords="%(x1)d,%(y1)d,%(x2)d,%(y2)d" '
                             'alt="%(title)s" title="%(title)s" href="%(href)s"%(target_attr)s />\n' % vals)
                f.write('</map>\n')
            f.write('<iframe name="empty" height="0" width="0" src="about:blank"></iframe>\n')
            f.write('<div><code>' + str(self.args) + '</code></div></body></html>\n')


def main():
    """Program entry point."""
    parser = argparse.ArgumentParser(description='Plot ROI file.')
    ngs_roi.argparse.addFileArguments(parser)
    ngs_roi.argparse.addPlotGridArguments(parser)
    ngs_roi.argparse.addLinkArguments(parser)
    args = parser.parse_args()
    ngs_roi.argparse.applyFileDefaults(args)

    app = PlotThumbnailsApp(args)
    return app.run()


if __name__ == '__main__':
    sys.exit(main())
