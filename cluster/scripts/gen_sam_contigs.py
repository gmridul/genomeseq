#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import getopt
import re
import math
import time

regExStr = '([0-9]+)([MIDNSHPX=])'
regExS = '(.*?[0-9]+)([ACGRTY]+)'
reExM = '[0-9]+'

def usage():
    print '%s -f /path/to/sam-file' % (sys.argv[0])


def processLine(line):
    elts = line.split()
    startPos = currentPos = int(elts[3])
    cigarString = elts[5]
    if cigarString != '*':
        currentStr = cigarString
        reM = re.match(regExStr, currentStr)
        while reM != None:
            length = int(reM.group(1))
            mChar = reM.group(2)
            end = reM.end()
            if mChar != 'D':
                currentPos += length
            if end < len(currentStr):
                currentStr = currentStr[end:]
                reM = re.match(regExStr, currentStr)
            else:
                break
    return (startPos, currentPos-1)

def getFirstLine(f):
    for line in f:
        if line[0] == '@':
            continue
        return processLine(line)

def process(fileName):
    with open(fileName, 'r') as f:
        (currentStart, currentEnd) = getFirstLine(f)
        lastStart = currentStart
        for line in f:
            (start, end) = processLine(line)
            #print "[", start, end, "]"
            assert(lastStart <= start)
            lastStart = start
            if (start > currentEnd+1):
                print currentStart, currentEnd
                currentStart = start
            if (end > currentEnd):
                currentEnd = end
        print currentStart, " ", currentEnd
        


def main(argv):
    filename = ""
    try:
      (opts, args) = getopt.getopt(argv, 'f:h', ['file=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for (opt, arg) in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-f', '--file'):
            filename = arg
        else:
            usage()
            sys.exit()
    if (filename == ""):
        usage()
        sys.exit(-1)
    process(filename)

if __name__ == '__main__':
    main(sys.argv[1:])
