#!/usr/bin/env python2

class RoiRecord(object):
    """Represent one record in a ROI file."""

    def __init__(self, ref, start_pos, end_pos, region_name, region_length,
                 strand, max_count, data, points):
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

    def __str__(self):
        return 'RoiRecord(%s, %s, %s, %s, %s, %s, %s, %s, len([...])==%s)' % \
                   (repr(self.ref), self.start_pos, self.end_pos,
                    self.region_name, self.region_length, repr(self.strand),
                    self.max_count, self.data, len(self.points))
    def __repr__(self):
        return self.__str__()


def loadRoi(path, max_count=0):
    """Load ROI file and return it as a list of RoiRecord objects.

    NA values are translated to 0.
    """
    data_keys = []
    result = []
    i = 0
    with open(path, 'rb') as f:
        for line in f:
            if line.startswith('##'):
                data_keys = line[2:].split('\t')[7:-1]
            if line[0] == '#':
                continue
            if max_count > 0 and i >= max_count:
                break
            i += 1
            vals = line.split()
            region_length = int(vals[4])
            data = vals[7:-1]
            points = [int(x) for x in vals[-1].split(',')]
            r = RoiRecord(vals[0], int(vals[1]) - 1, int(vals[2]), vals[3],
                         region_length, vals[5], int(vals[6]), data, points)
            result.append(r)
    #print '  => Loaded %d records.' % len(result)
    return data_keys, result
