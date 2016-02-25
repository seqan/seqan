#!/usr/bin/env python2
"""Input / Output Routines for ROI files.

You can load the records sequentially one by one by iterating over a
RoiFile object:

f1 = ngs_roi.io.RoiFile('in.roi')
for r in f1:
    print r.ref, r.start_pos, r.end_pos

f2 = ngs_roi.io.RoiFile('in.roi')
print f1.next().ref
print f1.next().ref
print f1.next().ref

Alternatively, you can load all or a subset of the records:

ten_records = ngs_roi.io.load('in.roi', 10)  # first ten
all_records = ngs_roi.io.load('in.roi')
"""

__author__    = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'
__copyright__ = 'Copyring 2013, Freie Universitaet Berlin'
__license__   = 'BSD 3-clause'


class RoiRecord(object):
    """Represent one record in a ROI file.

    :ivar ref:name of reference
    :type ref:str
    :ivar start_pos:0-based start position
    :type start_pos:int
    :ivar end_pos:0-based end position
    :type end_pos:int
    :ivar region_name:name of the region
    :type region_name:str
    :ivar region_length:length of the region
    :type region_length:int
    :ivar strand:strand of the region, one of ['+', '-']
    :type strand:str
    :ivar max_count:highest coverage
    :type max_count:int
    :ivar data:values of the extended data fields
    :type data:list of str
    :ivar points:the coverage over the region length
    :type points:list of ints

    :ivar data_keys:list with key names for the data fields
    :type data_keys:list of str or None
    """

    def __init__(self, ref, start_pos, end_pos, region_name, region_length,
                 strand, max_count, data, points, data_keys=None):
        """Initialize RoiRecord."""
        self.ref = ref
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.strand = strand
        self.region_length = region_length
        self.region_name = region_name
        self.max_count = max_count
        self.data = data
        self.points = points
        self.data_keys = data_keys

    def __str__(self):
        return 'RoiRecord(%s, %s, %s, %s, %s, %s, %s, %s, len([...])==%s, %s)' % \
                   (repr(self.ref), self.start_pos, self.end_pos,
                    self.region_name, self.region_length, repr(self.strand),
                    self.max_count, self.data, len(self.points), self.data_keys)
    def __repr__(self):
        return self.__str__()


class RoiFile(object):
    """File of ROI records.

    Can be used as an iterator.

    :ivar data_keys:Keys of additional data
    :type data_keys:list of str
    """

    def __init__(self, path):
        """Open ROI file and load header."""
        self.path = path
        self.f = open(path, 'rb')
        self.line = None
        self.data_keys = self._loadKeys(path)

    def _loadKeys(self, path):
        """Load keys from header, skipping comments."""
        while True:
            self.line = self.f.readline()
            if self.line == '':  # EOF
                return None
            if self.line.startswith('##'):
                keys = self.line[1:].strip().split('\t')
                self.line = self.f.readline()
                return keys[7:-1]
            if self.line.startswith('#'):
                continue  # comment
            break

    def __iter__(self):
        """Return iterator (self)."""
        return self

    def next(self):
        """Return next record."""
        if self.line == '':  # EOF
            raise StopIteration
        l = self.line
        self.line = self.f.readline()
        return self._buildRecord(l)

    def _buildRecord(self, line):
        """Build RoiRecord from ROI line."""
        vals = line.split()
        region_length = int(vals[4])
        data = vals[7:-1]
        points = [int(x) for x in vals[-1].split(',')]
        return RoiRecord(vals[0], int(vals[1]) - 1, int(vals[2]), vals[3],
                         region_length, vals[5], int(vals[6]), data, points,
                         self.data_keys)



def load(path, max_count=0):
    """Load ROI file and return it as a list of RoiRecord objects.

    NA values are translated to 0.
    """
    result = []
    for i, x in enumerate(RoiFile(path)):
        if max_count > 0 and i >= max_count:
            break
        result.append(x)
    return result
