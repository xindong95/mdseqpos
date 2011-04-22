#!/usr/bin/env python
# weili@jimmy.harvard.edu


import os.path, sys, re
from numpy import *

'''
All fields from selected table in UCSC, the first line is the headings.
The program will replace all "." to "_" and remove all "#" in the headings.
http://genome.ucsc.edu/cgi-bin/hgTables?org=Human&db=hg17&hgsid=63296757

float numpy array for all numerical entry
regular python list for all exon.*s entry
string array for all the remaining
'''

class UCSCTable(object):
    ''' UCSCTable  object.
    UCSCTable(filename,  *pars): like "chrom", "chromStart", if pars are not specified,
    it will read everything.
    Return: self.head = list
    '''
    def __init__(self, filename,  *pars):
        print >> sys.stderr, 'Loading annotations'
        
        self.filename = os.path.basename(filename)
        file = open(filename)
        headline = re.sub('\.', '_', file.readline().strip())
        headline = re.sub('#', '', headline)
        headline = headline.split('\t')
        index = []

        if len(pars) == 0 :
            self.head = headline
            index = range(len(self.head))
        else:
            self.head = []
            for i, x in enumerate(headline):
                if x in pars:
                    self.head.append(x)
                    index.append(i)
        self.head = char.array(self.head)

        for x in  self.head:
            self.__setattr__(x, [])

        for iline, line in enumerate(file):
            line = line.strip().split('\t')
            for i, x in enumerate(self.head):
                self.__getattribute__(x).append(line[index[i]])

        for x in self.head:
            try:
                int(self.__getattribute__(x)[0])
                self.__setattr__(x, array([int(i) for i in self.__getattribute__(x)],dtype=int))
            except:
                if not re.search("exon.*s", x):
                    self.__setattr__(x,char.array(self.__getattribute__(x)))
                    self.__setattr__(x,self.__getattribute__(x))

