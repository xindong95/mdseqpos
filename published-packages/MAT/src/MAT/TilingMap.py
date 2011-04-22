#!/usr/bin/env python
# weili@jimmy.harvard.edu

'''
Maping annotation from UCSC Genome Browser to bed regions
For Gene-like annotation, it should have 'strand', 'txStart', 'txEnd', 'exonStarts', 'exonEnds'
For non-gene annotation , it should have chromStart, chromEnd

Now support:

AceViewGenes
EnsemblGenes
KnownGenes
RefSeqGenes
NHGRI_DNasel-HS
sno_miRNA
'''

import re, sys, UCSCTable, os.path, time
import cPickle as pickle
from numpy import *


class TilingMap(object):
    '''TilingMap(UCSCTable.filename):
    '''
    GeneAnno = ['AceViewGenes', 'EnsemblGenes', 'KnownGenes', 'RefSeqGenes']
    NonGeneAnno = ['NHGRI_DNasel-HS', 'sno_miRNA']

    def __init__(self, filename):
        self.__dict__ = UCSCTable.UCSCTable(filename).__dict__
        headLength = self.head.shape[0]
        print >> sys.stderr, 'Starting annotation'
        
        # normal entries
        for x in ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'exonStarts', 'exonEnds',\
                  'chromStart', 'chromEnd']:
            hit = [i for i in range(headLength) if re.search(x+'$',self.head[i])]
            if not len(hit):
                continue
            fullname = self.head[hit[0]]
            self.__setattr__(x, self.__getattribute__(fullname))
            if fullname != x:
                self.__delattr__(fullname)

        # refSeq
        if self.filename == 'RefSeqGenes':
            refLink_name = self.head[[i for i in range(headLength) if re.search('refLink_name',self.head[i])][0]]
            refLink_product = self.head[[i for i in range(headLength) if re.search('refLink_product',self.head[i])][0]]
            self.name = (self.name + ': ' +self.__getattribute__(refLink_name)).strip() \
            + (',' + self.__getattribute__(refLink_product)).strip()

        # knownGene
        if self.filename == 'KnownGenes':
            kgXref_geneSymbol = self.head[[i for i in range(headLength) if re.search('kgXref_geneSymbol',self.head[i])][0]]
            kgXref_description = self.head[[i for i in range(headLength) if re.search('kgXref_description',self.head[i])][0]]
            self.name = (self.name + ': ' +self.__getattribute__(kgXref_geneSymbol)).strip() + \
            (',' + self.__getattribute__(kgXref_description)).strip()

        # EnsemblGenes
        if self.filename == 'EnsemblGenes':
            sfDescription_description = self.head[[i for i in range(headLength) if re.search('sfDescription_description',self.head[i])][0]]
            self.name = self.name + ': ' +self.__getattribute__(sfDescription_description).strip()

        # sno_miRNA
        if self.filename == 'sno_miRNA':
            type = self.head[[i for i in range(headLength) if re.search('type',self.head[i])][0]]
            self.name = self.name + ': ' +self.__getattribute__(type).strip()

    def __sc(self, n) :
        if n < 0 :
            n = -1 * n
            head = '-'
        else:
            head = ''

        if n <1000: return head+str(int(n))
        else:
            if n%1000 < 10 : return head+"%d,00%d" % (n/1000, n%1000)
            elif n%1000 <100 : return head+"%d,0%d" % (n/1000, n%1000)
            else : return head+ "%d,%d" % (n/1000, n%1000)


    def MapFromSites(self, bedfile, out = ''):
        '''MapFromSites(bedfile, out = ''): Mapping bed regions to non-gene-like annotation
        '''
        stdout = sys.stdout
        if out:
            sys.stdout = open(out, 'w', 0)    # out file
        print 'Distance\tName\t'

        oldchr = 999
        for line in open(bedfile):
            if line[:3] != 'chr':
                continue
            x = line.split()
            chr, pos = x[0], (int(x[1]) + int(x[2]))/2.0
            if oldchr != chr:
                segment = where(self.chrom == chr)[0]
                chromPos = (self.chromStart[segment] + self.chromEnd[segment]) / 2.0
                oldchr = chr

            if not segment.shape[0]:             # no annotation on this chromosome
                continue

            upstream = abs(chromPos - pos).argmin()
            upstream_dis = abs(chromPos - pos)[upstream]
            upstream = segment[upstream]
            print "%s\t%s\t" % (self.__sc(upstream_dis), self.name[upstream]),
            print line.strip()

        sys.stdout = stdout


    def MapGeneFromSites(self, bedfile, out = ''):
        '''MapGeneFromSites(self, bedfile, out = ''): Mapping bed regions to gene-like annotation
        '''

        stdout = sys.stdout
        if out:
            sys.stdout = open(out, 'w', 0)    # out file
        print '5Prime_dis\t5Prime_Gene\t3Prime_dis\t3Prime_Gene\tMiddle_Anno\tMiddle'

        oldchr = 999
        for line in open(bedfile):
            if line[:3] != 'chr':
                continue
            x = line.split()
            chr, pos = x[0], (int(x[1]) + int(x[2]))/2          # pos: middle if region
            if oldchr != chr:
                segment = where(self.chrom == chr)[0]
                strand = self.strand[segment]
                name = self.name[segment]
                txStart = self.txStart[segment]
                txEnd = self.txEnd[segment]
                #txStart = self.txStart[segment]*(strand == '+') + self.txEnd[segment]*(strand == '-')
                #txEnd = self.txStart[segment]*(strand == '-') + self.txEnd[segment]*(strand == '+')
                oldchr = chr
                
            if not segment.shape[0]:             # no annotation on this chromosome
                print " \t" * 6, line.strip()
                continue

            # to the closest 5' or 3'.
            for tx in [txStart*(strand == '+') + txEnd*(strand == '-'), \
                       txStart*(strand == '-') + txEnd*(strand == '+')]:
                dis = tx - pos
                min_id = abs(dis).argmin()
                if (dis[min_id] <= 0 and strand[min_id] == '-') or\
                    (dis[min_id] >=0 and strand[min_id] == '+'):
                    dis_min = abs(dis)[min_id]
                else:
                    dis_min = -1 * abs(dis)[min_id]
                print  self.__sc(dis_min), '\t', name[min_id],'\t',

            # middle
            mid_id = where((pos > txStart) * (pos < txEnd) )[0]
            if not mid_id.shape[0]:                   # not in the middle of annotation
                print  ' \t \t',
            else:
                mid_id = segment[mid_id[0]]         # the first gene
                exonEnds = array([int(x) for x in self.exonEnds[mid_id][:-1].split(',')],dtype=int)
                exonStarts = array([int(x) for x in self.exonStarts[mid_id][:-1].split(',')],dtype=int)
                exonEnds = where(exonEnds >= pos)[0]
                exonStarts = where(exonStarts >= pos)[0]
                if not exonStarts.shape[0] or exonStarts[0] > exonEnds[0]:
                    print "E%d\t%s\t" % (exonEnds[0]+1, self.name[mid_id]),
                else:
                    print "I%d\t%s\t" % (exonEnds[0], self.name[mid_id]),

            print line.strip()
        sys.stdout = stdout


    def MapGeneFromSites_Old(self, bedfile, out = ''):
        '''MapGeneFromSites(self, bedfile, out = ''): Mapping bed regions to gene-like annotation
        '''

        stdout = sys.stdout
        if out:
            sys.stdout = open(out, 'w', 0)    # out file
        print 'Upstream_dis\tUpstream\tDownstream_dis\tDownstream\tMiddle_Anno\tMiddle'

        oldchr = 999
        for line in open(bedfile):
            if line[:3] != 'chr':
                continue
            x = line.split()
            chr, pos = x[0], (int(x[1]) + int(x[2]))/2          # pos: middle if region
            bw = (int(x[2]) - int(x[1]))/2                      # bandwidth of the region
            if oldchr != chr:
                segment = where(self.chrom == chr)[0]
                txEnd = self.txEnd[segment]
                txStart = self.txStart[segment]
                strand = self.strand[segment]
                oldchr = chr

            if not segment.shape[0]:             # no annotation on this chromosome
                print " \t" * 6, line.strip()
                continue

            # upstream, bed region in the leftplus of annotation, etc.
            leftplus = where((txStart >= (pos-bw)) *  (strand == '+'))[0]
            rightminus = where((txEnd <= (pos+bw)) *  (strand == '-'))[0]
            # downstram
            leftminus = where((txStart >= (pos-bw)) *  (strand == '-'))[0]
            rightplus = where((txEnd <= (pos+bw)) *  (strand == '+'))[0]
            for (left, right) in [(leftplus, rightminus), (leftminus, rightplus)]:
                if left.shape[0] and right.shape[0]:          # find two genes
                    left = left[0]
                    right = right[-1]
                    if (txStart[left] - pos) < (pos - txEnd[right]):
                        upstream = segment[left]                # id for the nearest gene
                        upstream_dis = txStart[left] - pos
                    else:
                        upstream = segment[right]
                        upstream_dis = pos - txEnd[right]
                elif left.shape[0] :                            # only left
                    left = left[0]
                    upstream = segment[left]
                    upstream_dis = txStart[left] - pos
                elif right.shape[0]:                           # only right
                    right = right[-1]
                    upstream = segment[right]
                    upstream_dis = pos - txEnd[right]
                else:
                    print " \t \t",
                    continue
                print  self.__sc(upstream_dis), '\t', self.name[upstream],'\t',

            # middle
            upstream = where(((pos-bw) > txStart) * ((pos+bw) < txEnd) )[0]
            if not upstream.shape[0]:                   # not in the middle of annotation
                print  ' \t \t',
            else:
                upstream = segment[upstream[0]]         # the first gene
                exonEnds = array([int(x) for x in self.exonEnds[upstream][:-1].split(',')],dtype=int)
                exonStarts = array([int(x) for x in self.exonStarts[upstream][:-1].split(',')],dtype=int)
                exonEnds = where(exonEnds >= pos)[0]
                exonStarts = where(exonStarts >= pos)[0]
                if not exonStarts.shape[0] or exonStarts[0] > exonEnds[0]:
                    print "E%d\t%s\t" % (exonEnds[0]+1, self.name[upstream]),
                else:
                    print "I%d\t%s\t" % (exonEnds[0], self.name[upstream]),

            print line.strip()
        sys.stdout = stdout



    def Mapping(self, bedfile, out = ''):
        '''Mapping(self, bedfile, out = ''):
        '''
        print >> sys.stderr, 'Mapping Annotation ', self.filename, time.asctime()

        if self.filename in self.GeneAnno :
            self.MapGeneFromSites(bedfile, out)
        elif self.filename in self.NonGeneAnno:
            self.MapFromSites(bedfile, out)



#read program options
if __name__ == '__main__':
    if len(sys.argv) <2:
        print "Usage: %prog annotation_file  bedfile"
        sys.exit(0)

    ref = TilingMap(sys.argv[1])
    for bedfile in sys.argv[2:]:
        ref.Mapping(bedfile, bedfile+ '.' + os.path.basename(sys.argv[1]) + '.xls')
