#!/usr/bin/env python2
"""Support code for writing the NGS ROI report generation apps."""

__author__    = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'
__copyright__ = 'Copyring 2013, Freie Universitaet Berlin'
__license__   = 'BSD 3-clause'


import os
import os.path
import subprocess


class App(object):
    """Base class with helper functionality for NGS ROI report generators."""

    def __init__(self, args):
        self.args = args

    def prepareOutDir(self):
        """Create output directory if necessary."""
        if not os.path.exists(self.args.out_dir):
            os.makedirs(self.args.out_dir)

    def buildHref(self, ref, begin_pos, end_pos):
        vals = {'ref': ref, 'start_pos': begin_pos + 1, 'end_pos': end_pos}
        if self.args.link_type == 'local_igv':
            vals['host'] = self.args.igv_host
            vals['port'] = self.args.igv_port
            return 'http://%(host)s:%(port)d/goto?locus=%(ref)s:%(start_pos)d-%(end_pos)d' % vals
        else:  # self.args.link_type == 'ucsc'
            vals['org'] = self.args.ucsc_org
            vals['db'] = self.args.ucsc_db
            vals['ref'] = self.args.ucsc_chr_prefix + vals['ref']
            return 'http://genome.ucsc.edu/cgi-bin/hgTracks?org=%(org)s&db=%(db)s&position=%(ref)s:%(start_pos)d-%(end_pos)d' % vals


class PlotThumbnailsRunner(object):
    """Helper class for running the roi_plot_thumbnails program.

    :ivar args: Arguments as returned by argparse.
    """

    def __init__(self, args):
        self.args = args

    def run(self):
        cmd_args = ['-if', self.args.in_file,
                    '--input-file-file-ext', 'roi',
                    '-o', os.path.join(self.args.out_dir, 'thumbnail_'),
                    '--max-rois', self.args.max_rois,
                    '--max-value', self.args.max_value,
                    '--num-cols', self.args.num_cols,
                    '--num-rows', self.args.num_rows,
                    '--plot-height', self.args.plot_height,
                    '--plot-width', self.args.plot_width,
                    '--border-width', self.args.border_width,
                    '--spacing', self.args.spacing]
        cmd_args = ['roi_plot_thumbnails'] + map(str, cmd_args)
        #import pdb; pdb.set_trace()
        import sys
        print >>sys.stderr, 'Running %s' % ' '.join(cmd_args)
        p = subprocess.Popen(cmd_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        res = p.wait()
        if res:
            print 'ERROR', p.stdin, p.stderr
        return res
