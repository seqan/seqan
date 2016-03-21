#!/usr/bin/env python2
"""Helper for ROI tools when using argparse module.

This module contains helper functions for setup of argparse.ArgumentParser
that is common to multiple report generating apps.
"""

__author__    = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'
__copyright__ = 'Copyring 2013, Freie Universitaet Berlin'
__license__   = 'BSD 3-clause'


import os.path


def addFileArguments(parser):
    """Adds --in-file, --out-file, and --out-dir to argument parser.

    These parameters are used for the input ROI file, the output HTML file and
    the directory to write supplementary files to.
    """
    group = parser.add_argument_group('Input / Output',
                                      'Input and output related parameters.')
    group.add_argument('--in-file', metavar='ROI', dest='in_file',
                       required=True, help='ROI file to read')
    group.add_argument('--out-file', dest='out_file', required=True,
                       help='path to output HTML file')
    group.add_argument('--out-dir', dest='out_dir',
                       help='directory to write supplementary files to; '
                       'defaults to path of --out-file.')


def applyFileDefaults(args):
    """Apply file-related default values (copying paths if not set)."""
    if not args.out_dir:
        args.out_dir = os.path.dirname(args.out_file) or '.'


def addPlotGridArguments(parser, default_plot_height=30, default_plot_width=30):
    """Adds arguments related to plot grids.

    This is used for the thumbnail plotting but also for the tables.
    """
    group = parser.add_argument_group('Plot Grid Configuration',
                                      'Arguments for the plot image grid.')

    group.add_argument('--max-rois', dest='max_rois', metavar='NUM',
                        type=int, default=0,
                        help='Maximal number of ROIs, 0 for all.')
    group.add_argument('--max-value', dest='max_value', metavar='NUM',
                        type=int, default=0,
                        help='Largest y value to plot, 0 for all.')

    group.add_argument('--num-rows', dest='num_rows', metavar='ROWS',
                        type=int, default=50,
                        help='Number of rows per grid.')
    group.add_argument('--num-cols', dest='num_cols', metavar='COLS',
                       type=int, default=40,
                       help='Number of columns per grid.')
    
    group.add_argument('--plot-height', dest='plot_height', metavar='HEIGHT',
                       type=int, default=default_plot_height, help='Height of one plot in px.')
    group.add_argument('--plot-width', dest='plot_width', metavar='WIDTH',
                       type=int, default=default_plot_width, help='Width of one plot in px.')

    group.add_argument('--border-width', dest='border_width', metavar='WIDTH',
                       type=int, default=0, help='Border width.')
    group.add_argument('--spacing', dest='spacing', metavar='SPACING',
                       type=int, default=2, help='Spacing.')


def addLinkArguments(parser):
    """Adds arguments related to link creation.

    These parameters control the created links to local IGV browser HTTP
    remote control or the UCSC genome browser.
    """
    group = parser.add_argument_group('HTML Links', 'Arguments for HTML link creation.')

    group.add_argument('--link-target', dest='link_target', metavar='TARGET',
                       default='_blank', choices=['_blank', '_top'],
                       help='Select the link target to create (_blank or _top).')
    
    group.add_argument('--link-type', dest='link_type', metavar='TARGET',
                       default='local_igv', choices=['local_igv', 'ucsc'],
                       help='Select the type of links to create.  One of '
                       '"local_igv" and "ucsc".')
    
    group.add_argument('--igv-host', dest='igv_host', metavar='HOST',
                       default='localhost', help='Host for IGV link.')
    group.add_argument('--igv-port', dest='igv_port', metavar='PORT',
                       type=int, default='60151', help='Port for IGV link.')

    group.add_argument('--ucsc-org', dest='ucsc_org', metavar='ORG',
                       default='human', help='Organism for UCSC browser link.')
    group.add_argument('--ucsc-db', dest='ucsc_db', metavar='DB',
                       default='hg18', help='Assembly version for UCSC browser link.')
    group.add_argument('--ucsc-chr-prefix', dest='ucsc_chr_prefix', metavar='PREFIX',
                       default='', help='Prefix for chromosome names in UCSC browser.')
