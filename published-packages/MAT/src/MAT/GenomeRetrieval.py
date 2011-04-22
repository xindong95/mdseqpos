#!/usr/bin/env python
# weili@jimmy.harvard.edu

import sys, os.path, re, time


class GenomeRetrieval(object):
    def __init__(self, genomedir = \
        "/misc/iris/acct/weili/database/humanhg17_May2004/masked/*.fa.masked.simple" ):
        self.genomedir = genomedir


    def GetChrSeq(self, chr):
        seq=[]
        for x in open(os.path.dirname(self.genomedir) + os.path.sep + re.sub(\
                    '\*', chr, os.path.basename(self.genomedir))):
            if x[0] == '>':
                continue
            seq.append(x.strip())
        return ''.join(seq)


    def GetBedSeq(self, bedname, seqname = ''):
        print >> sys.stderr, 'Retrieval sequences ... ', time.asctime()
        oldchr = ''
        if seqname:
            output = open(seqname, 'w', 0)
        else:
            output = sys.stdout

        for line in open(bedname) :
            if line[:3] != 'chr':
               continue
            x=line.split()
            chr, left, right = x[0], int(x[1]), int(x[2])
            print >>output, ">%s" % (line),
            if oldchr != chr:                   # new chr
                seq = self.GetChrSeq(chr)
                oldchr = chr
            print >>output, seq[left-1:right]


#read program options
if __name__ == '__main__':
    if len(sys.argv) <2:
        print "Usage: %prog bedfile"
        sys.exit(0)

    genome = GenomeRetrieval()
    genome.GetBedSeq(sys.argv[1])




