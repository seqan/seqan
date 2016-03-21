#!/usr/bin/env python2
"""Convert SeqAn profiling information into PDF graphic.

USAGE: profile2pdf.py <program.profile.txt> <out.pdf>
"""

from __future__ import with_statement

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import math
import sys

import cairo


IGNORE_LIMIT = 0.00001

# Automatic legend colors.
COLORS = [
  "#3366FF",
  "#FF33CC",
  "#FF6633",
  "#CCFF33",
  "#33FF66",
  "#33CCFF",
  "#002EB8",
  "#B88A00",
  "#CC33FF",
  "#FF3366",
  "#FFCC33",
  "#66FF33",
  "#33FFCC",
  "#003DF5",
  "#F5B800"
  ]
def htmlColorToRgb(colorstring):
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError, "input #%s is not in #RRGGBB format" % colorstring
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]
    return (r / 255.0, g / 255.0, b / 255.0)
COLORS = map(htmlColorToRgb, COLORS)  

class Meta(object):
  def __init__(self, beginTimestamp, endTimestamp):
    self.beginTimestamp = beginTimestamp
    self.endTimestamp = endTimestamp

class JobType(object):
  """Describe a job type."""
  def __init__(self, identifier, shortName, longName=None, color=None):
    self.identifier = identifier
    self.shortName = shortName
    self.longName = longName or shortName
    self.color = color or COLORS[identifier % len(COLORS)]
  
  @classmethod
  def fromString(klass, s):
    columns = s.split('\t')
    if columns[0] != '@EVENT':
      print >>sys.stderr, 'First column\'s value was not "@EVENT@'
      sys.exit(1)
    identifier = int(columns[1])
    shortName = columns[2]
    longName = columns[3]

    # Read in optional arguments.
    color = None
    for col in columns[4:]:
      key, value = col.split(':')
      if key == 'COLOR':
        color = value
    return JobType(identifier, shortName, longName=longName, color=color)

class Event(object):
  """Describes an event."""
  def __init__(self, threadId, isBegin, jobType, timestamp):
    self.threadId = threadId
    self.isBegin = isBegin
    self.jobType = jobType
    self.timestamp = timestamp

  @classmethod
  def fromString(klass, s):
    columns = s.split('\t')
    threadId = int(columns[0])
    if columns[1] not in ['BEGIN', 'END']:
      print >>sys.stderr, 'Second column\'s value was not BEGIN or END'
      sys.exit(1)
    isBegin = columns[1] == 'BEGIN'
    jobType = int(columns[2])
    timestamp = float(columns[3])
    return Event(threadId, isBegin, jobType, timestamp)


class Section(object):
  """Describe a section in the program run."""
  def __init__(self, jobType, beginTime, endTime, parent=None):
    self.children = []
    self.jobType = jobType
    self.beginTime = beginTime
    self.endTime = endTime
    self.parent = parent
    if self.parent:
      self.parent.children.append(self)

def buildSections(events):
  forest = []
  sections = []
  stack = []
  for e in events:
    if e.isBegin:
      if not stack:
        section = Section(e.jobType, e.timestamp, endTime=None, parent=None)
        forest.append(section)
      else:
        section = Section(e.jobType, e.timestamp, endTime=None, parent=stack[-1])
      sections.append(section)
      stack.append(section)
    else:  # e.isBegin
      assert stack
      assert stack[-1].jobType == e.jobType
      section = stack[-1]
      section.endTime = e.timestamp
      stack.pop()
  return forest, sections

def printSection(section, jobTypes, offset, level=0):
  span = section.endTime - section.beginTime
  print '%s%s %f  (%f to %f)' % ('\t' * level, jobTypes[section.jobType].shortName, span, section.beginTime - offset, section.endTime - offset)
  for s in section.children:
    printSection(s, jobTypes, offset, level+1)

def loadFile(path):
  with open(path, 'r') as f:
    line = f.readline()
    if line.strip() != '@SQN:PROFILE':
      print >>sys.stderr, 'Invalid file, does not start with "@SQN:PROFILE"'
      sys.exit(1)
    line = f.readline()
    if not line.startswith('@TIME'):
      print >>sys.stderr, 'Invalid file, second line does not start with "@TIME"'
      sys.exit(1)
    meta = Meta(*map(float, line.strip().split('\t')[1:]))
    # Load job types.
    jobTypes = []
    while True:
      line = f.readline()
      if not line or not line.startswith('@EVENT'):
        break  # End of file, no more job types.
      jobTypes.append(JobType.fromString(line.strip()))
    # Events.
    events = []
    while True:
      if not line:
        break
      events.append(Event.fromString(line.strip()))
      line = f.readline()
    # Remove redundant entries.
    if False:
      events2 = []
      for i in range(len(events)):
        if not events2 or not i + 1 < len(events):
          events2.append(events[i])
          continue
        while True:
          if not events2[-1].isBegin and events[i].isBegin and events[-1].jobType == events[i].jobType and events[-1].threadId == events[i].threadId and events[-1].timestamp - events[i].timestamp < 0.0001 and not events[i + 1].isBegin and events[i].jobType == events[i + 1].jobType and events[i].threadId == events[i + 1].threadId:
            i += 2
          else:
            break
        events2.append(events[i])
  return meta, jobTypes, events

POINTS_SPACE_OUTER = 10
POINTS_PER_SECOND = 10
POINTS_BAR_HEIGHT = 5
POINTS_SPACE = 2
POINTS_KEY_ENTRY_HEIGHT = 5;

def drawBox(cr, jobTypes, section, offset, threadId, level):
  assert level < 10
  x = POINTS_SPACE_OUTER + (section.beginTime - offset) * POINTS_PER_SECOND
  y = POINTS_SPACE_OUTER + POINTS_SPACE * threadId + POINTS_BAR_HEIGHT * threadId + level * 0.1 * POINTS_BAR_HEIGHT
  width = (section.endTime - section.beginTime) * POINTS_PER_SECOND
  height = (1.0 - 0.1 * level) * POINTS_BAR_HEIGHT
  #print 'rectangle(%s, %s, %s, %s), level = %s' % (x, y, width, height, level)
  cr.set_source_rgb(*jobTypes[section.jobType].color)
  cr.rectangle(x, y, width, height)
  cr.fill()
  cr.set_line_width(0.2)
  cr.set_source_rgb(0, 0, 0)
  cr.move_to(x + width, y)
  cr.line_to(x + width, y + height)
  cr.stroke()

def drawBoxesForSection(cr, jobTypes, section, offset, threadId, level=0):
  drawBox(cr, jobTypes, section, offset, threadId, level)
  for s in section.children:
    drawBoxesForSection(cr, jobTypes, s, offset, threadId, level + 1)

def drawKey(cr, jobTypes, threadCount):
  for i, jobType in enumerate(jobTypes):
    x = POINTS_SPACE_OUTER
    y = POINTS_BAR_HEIGHT * threadCount + POINTS_SPACE_OUTER + POINTS_SPACE * (threadCount + 1) + POINTS_KEY_ENTRY_HEIGHT * i
    width = POINTS_KEY_ENTRY_HEIGHT * 2
    height = POINTS_KEY_ENTRY_HEIGHT
    cr.set_source_rgb(*jobTypes[i].color)
    cr.rectangle(x, y, width, height)
    cr.fill()
    cr.set_source_rgb(0, 0, 0)
    cr.set_font_size(POINTS_KEY_ENTRY_HEIGHT)
    cr.move_to(x + 3 * POINTS_KEY_ENTRY_HEIGHT, y + 0.8 * POINTS_KEY_ENTRY_HEIGHT)
    cr.show_text(jobType.shortName)

def drawScale(cr, totalBegin, totalEnd, threadCount):
  cr.set_line_width(0.2)
  cr.set_font_size(POINTS_KEY_ENTRY_HEIGHT * 0.5)
  for i in range(0, int(totalEnd-totalBegin) + 1):
    # Draw ticks at top.
    cr.set_source_rgb(0, 0, 0)
    cr.move_to(POINTS_SPACE_OUTER + POINTS_PER_SECOND * i, POINTS_SPACE_OUTER)
    if i % 5:  # small tick
      cr.line_to(POINTS_SPACE_OUTER + POINTS_PER_SECOND * i, 0.9 * POINTS_SPACE_OUTER);
    else:  # large tick
      cr.line_to(POINTS_SPACE_OUTER + POINTS_PER_SECOND * i, 0.8 * POINTS_SPACE_OUTER);
    cr.stroke()
    # Draw grid.
    cr.set_source_rgba(.3, .3, .3, 0.5)
    cr.move_to(POINTS_SPACE_OUTER + POINTS_PER_SECOND * i, POINTS_SPACE_OUTER)
    cr.line_to(POINTS_SPACE_OUTER + POINTS_PER_SECOND * i, POINTS_SPACE_OUTER + POINTS_SPACE * (threadCount - 1) + POINTS_BAR_HEIGHT * threadCount);
    cr.stroke()
    if not i % 5:  # at large ticks
      # Draw seconds display.
      cr.set_source_rgb(0, 0, 0)
      extents = cr.text_extents(str(i))
      cr.move_to(POINTS_SPACE_OUTER + POINTS_PER_SECOND * i - extents[2] / 2.0, 0.75 * POINTS_SPACE_OUTER);
      cr.show_text(str(i))

def breakDownTimesHelper(counter, section):
  counter[section.jobType] = counter.get(section.jobType, 0) + section.endTime - section.beginTime
  if section.parent:
    counter[section.parent.jobType] -= section.endTime - section.beginTime
  for child in section.children:
    breakDownTimesHelper(counter, child)

def breakDownTimes(jobTypes, forests):
  for threadId in sorted(forests.keys()):
    print 'Breakdown for thread #%d' % threadId
    counter = {}
    for section in forests[threadId]:
      breakDownTimesHelper(counter, section)
    for jobType in jobTypes:
      print '  %20s %10.5f' % (jobType.shortName, counter.get(jobType.identifier, 0))

def createDiagram(meta, jobTypes, forests, path):
  totalBegin = meta.beginTimestamp
  totalEnd = meta.endTimestamp
  totalTime = totalEnd - totalBegin
  # Create Cairo PDF surface.
  width = math.ceil(totalTime) * POINTS_PER_SECOND + POINTS_SPACE + 2 * POINTS_SPACE_OUTER
  height = POINTS_BAR_HEIGHT * len(forests) + POINTS_SPACE_OUTER + POINTS_SPACE * (len(forests) + 2) + POINTS_KEY_ENTRY_HEIGHT * len(jobTypes)

  cs = cairo.PDFSurface(path, width, height)
  cr = cairo.Context(cs)
  
  for threadId, forest in forests.iteritems():
    for section in forest:
      drawBoxesForSection(cr, jobTypes, section, totalBegin, threadId)
  drawKey(cr, jobTypes, len(forests))
  drawScale(cr, totalBegin, totalEnd, len(forests))
  cr.show_page()
  cs.finish()

def main(args):
  if len(args) != 3:
    print >>sys.stderr, 'Invalid number of arguments!'
    print >>sys.stderr, 'USAGE: profile2pdf.py <program.profile.txt> <out.pdf>'
    return 1
  
  # Load input file.
  print >>sys.stderr, 'Loading file', args[1]
  meta, jobTypes, events = loadFile(args[1])
  # Partition events by thread id.
  print >>sys.stderr, 'Partition events'
  eventsForThread = {}
  for e in events:
    eventsForThread.setdefault(e.threadId, []).append(e)

  # Build sections list and forest for each thread.
  print >>sys.stderr, 'Build sections'
  forests = {}
  sections = {}
  for threadId in sorted(eventsForThread.keys()):
    events = eventsForThread[threadId]
    f, s = buildSections(events)
    forests[threadId], sections[threadId] = f, s
    # Print sections (debug only):
    #print 'SECTIONS, threadId =', threadId
    #for x in f:
    #  printSection(x, jobTypes, s[0].beginTime)

  # Build diagram.
  print >>sys.stderr, 'Create diagram'
  createDiagram(meta, jobTypes, forests, args[2])

  # Show how much time each thread spent in each job type.
  breakDownTimes(jobTypes, forests)
  print 'TOTAL TIME: %f s' % (meta.endTimestamp - meta.beginTimestamp)
  
  return 0

if __name__ == '__main__':
  sys.exit(main(sys.argv))
