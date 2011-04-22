#!/usr/bin/env python
# weili@jimmy.harvard.edu

from UCSCTable import UCSCTable
import sys, time
from numpy import *
import cPickle as pickle

class Rep(object):
    '''Rep(filename = ''): load the repeatlib by filename or make a new one
    '''

    def __init__(self, filename = ''):
        self.repeatmsk = []
        self.simplerp = []
        self.segdup = []
        if  filename:
            self.__dict__ = pickle.load(open(filename))

        # to address the segmentation fault problem
#        tmp = self.segdup.__str__()

    def Make(self, *pars):
        '''Make(*pars): Make RepeatLib
        '''
        repeatmsk = UCSCTable(pars[0], 'genoName', 'genoStart', 'genoEnd')
        repeatmsk.chrom = repeatmsk.genoName
        simplerp = UCSCTable(pars[1], 'chrom', 'chromStart', 'chromEnd','copyNum')
        segdup = UCSCTable(pars[2], 'chrom', 'chromStart', 'chromEnd')

        #make ChrInd
        self.ChrInd  = []
        for x in [repeatmsk, simplerp, segdup]:
            for chr in x.chrom:
                if not chr in self.ChrInd:
                    self.ChrInd.append(chr)
        self.ChrInd.sort()

        # make matrix
        for rep in [repeatmsk, simplerp, segdup]:
            chrom_int = zeros(rep.chrom.shape,dtype=int8)
            for i, chr in enumerate(self.ChrInd):
                chrom_int[rep.chrom == chr] = i
            rep.chrom = chrom_int

        self.repeatmsk = transpose(array((repeatmsk.chrom, repeatmsk.genoStart, repeatmsk.genoEnd)))
        self.simplerp = transpose(array((simplerp.chrom, simplerp.chromStart, simplerp.chromEnd,
                                         simplerp.copyNum)))
        self.segdup = transpose(array((segdup.chrom, segdup.chromStart, segdup.chromEnd)))


    def MarkBed(self, bed, ChrInd, percent = 0.7):
        '''MarkBed( bed, ChrInd, percent = 0.7):
        bed: bed list
        percent: mark the region as repeat if over percentage of region is within the repeat lib
        '''
        print >> sys.stderr, 'Repeat Masking ... ', time.asctime()

        oldid = 999999999999
        for i, x in enumerate(bed):
            chr, st, end, name =  x[0],x [1], x[2], x[3]
            id = self.ChrInd.index(ChrInd[chr])
            # the new chromosome
            if id != oldid:
                repeatmsk = self.repeatmsk[self.repeatmsk[:,0] == id]
                simplerp = self.simplerp[self.simplerp[:, 0] == id]
                segdup = self.segdup[self.segdup[:, 0] == id]
                oldid = id

            # repeat hits
            repeatid = zeros(end-st,dtype=int)
            for rep in [repeatmsk, simplerp, segdup]:
                if rep.shape[0] == 0:   # no repeat in this chr
                    continue
                try:
                    repspan = ((rep[:,1] - end) * (rep[:,2] - st))
                    findrep = 0
                    if repspan.min() <= 0:   # find the repeat hit
                        for x, y in (rep[where(repspan <= 0)[0]][:, 1:3].astype(int)):
                            # big repeat unit
                            if (min(end, y) - max(x, st)) >= percent * (end - st):
                                findrep += 1
                            repeatid[arange(max(x, st) - st, min(end,y) -st)] = 1

                    if not findrep:
                        continue
                    if rep is simplerp:
                        bed[i][3] += '_Si'
                    if rep is segdup:
                        bed[i][3] += '_Se%d' % (findrep)
                except:
                    continue
            if (where(repeatid >0)[0].shape[0] * 1.0 / repeatid.shape[0] > percent):
                bed[i][3] = '_R' + bed[i][3]
        return bed



    def MarkBedFile(self, filename, percent = 0.7):
        '''MarkBed(self, bp, percent = 0.7):
        filename: bed file name
        percent: mark the region as repeat if the minmum percentage of region is within the repeat lib
        '''
        oldid = 999999999999
        for  line in open(filename):
            if line[:3] != 'chr':
                continue
            line = line.split()
            chr, st, end, name = line[0], int(line[1]), int(line[2]), line[3]
            id = self.ChrInd.index(chr)
            # the new chromosome
            if id != oldid:
                repeatmsk = self.repeatmsk[self.repeatmsk[:,0] == id]
                simplerp = self.simplerp[self.simplerp[:, 0] == id]
                segdup = self.segdup[self.segdup[:, 0] == id]
                oldid = id

            # repeat hits
            repeatid = zeros(end-st,dtype=int)
            for rep in [repeatmsk, simplerp, segdup]:
                if rep.shape[0] == 0:   # no repeat in this chr
                    continue
                repspan = ((rep[:,1] - end) * (rep[:,2] - st))
                findrep = 0
                if repspan.min() <= 0:   # find the repeat hit
                    for x, y in (rep[where(repspan <= 0)[0]][:, 1:3].astype(int)):
                        # big repeat unit
                        if (min(end, y) - max(x, st)) >= percent * (end - st):
                            findrep += 1
                        repeatid[arange(max(x, st) - st, min(end,y) -st)] = 1

                if not findrep:
                    continue
                if rep is simplerp:
                    name += '_Si'
                if rep is segdup:
                    name += '_Se%d' % (findrep)
            if (where(repeatid >0)[0].shape[0] * 1.0 / repeatid.shape[0] > percent):
                name += '_R' 
                line[3] = name
            print '\t'.join(line)


#read program options
# if __name__ == '__main__':
#     from optparse import OptionParser
#     import sys

#     parser = OptionParser(usage = "%prog [options] RepeatMsk simpleRepeat SegmentalDups", version = "%prog 02172008")

#     parser.add_option("-m",  action = "store_true", dest="make", default=False,
#                       help="make repeat lib")
#     parser.add_option("-b",  dest="bed",action = "store_true", default=False,
#                       help="mask bed files")

#     (options, args) = parser.parse_args()

#     if len(sys.argv) <2:
#         parser.print_help()
#         sys.exit()

#     if options.make:
#         rep = Rep()
#         rep.Make(args[0], args[1], args[2])
#         pickle.dump(rep.__dict__, open('RepeatLib','wb'), -1)

#     if options.bed:
#         rep = Rep('RepeatLib')
#         for file in args:
#             sys.stdout = open(file+'.rep', 'w', 0)
#             rep.MarkBedFile(file)

