#!/usr/bin/env python
# weili@jimmy.harvard.edu


from numpy import *
from Bed import Bed
from numpy import memmap as Map
import os.path, sys, cPickle
from glob import *

class phastCons(Bed):
    def __init__(self, bedfile, dir, win, ma, out, ymin = 0, ymax = 0):
        self.L = 17                 # number of letters per line
        Bed.__init__(self, bedfile)
        self.dir = dir
        self.win = win+ma
        self.ma = ma               # window for moving average

        try:
            from rpy import r
            r.pdf(out+'.ConservationPlot.pdf')
            self.Run()
            if not ymin:
                ymin = self.mscore[self.mscore>0].min()-0.05
            if not ymax:
                ymax = self.mscore.max()+0.05
            
            r.plot(range(-1*win, win+1), self.mscore, type = 'l', xlab = \
                   'Distance from the Center of Enriched Regions', \
                   ylab = 'Conservation Score', lwd= 3,  ylim = (ymin, ymax))
            r.dev_off()
        except:
            print >> sys.stderr, 'error import r using rpy, will not generate phastCons plot'
            print sys.exc_info()[0], sys.exc_info()[1]
    
            
    
    def Run(self):
        # init
        num = 0
        for chr in sorted(self.data.keys()):
            num += len(self.data[chr])
        self.score = zeros((num, (self.win)*2+1),dtype=float)

        num = 0
        offset = 0
        for chr in sorted(self.data.keys()):
            m = cPickle.load(open(self.dir + os.path.sep + chr + '.bin'))
            for line in self.data[chr]:
                print >> sys.stderr, line
                center = int((line[1] + line[0])/2)
                st, end = center - self.win, center + self.win +1
                for i, base in enumerate(xrange(st, end)):
                    try:                                    # out of range
                        self.score[num, i] = m[base]
                    except:
                        pass
                num += 1

        # moving average
        score = average(self.score, 0)
        self.mscore = zeros(score.shape,dtype=float)
        for x in range(self.ma, 2 * self.win + 1 - self.ma):
            self.mscore[x] = score[(x-self.ma): (x+self.ma+1)].mean()
        if self.ma :
            self.mscore = self.mscore[self.ma: self.ma*-1]/10000.0
        else:
            self.mscore = self.mscore / 10000.0
                
    
    
    def Old_Run(self):
        # init
        num = 0
        for chr in sorted(self.data.keys()):
            num += len(self.data[chr])
        self.score = zeros((num, (self.win)*2+1),dtype=float)

        num = 0
        offset = 0
        for chr in sorted(self.data.keys()):
            m = Map.open(self.dir + os.path.sep + chr, 'r')    # memory mapping object
            for line in self.data[chr]:
                print >> sys.stderr, line
                center = int((line[1] + line[0])/2)
                st, end = center - self.win, center + self.win +1
                for ind, base in enumerate(xrange(center, end)):
                    offset = m.find('\n'+str(base)+'\t', offset) + 1
                    if  offset :                            # found this one
                        break

                # get the entire region
                pos = []
                score = []
                for x in m[max(0,offset-self.L*(self.win+ind)) : \
                           offset+self.L*(self.win+1-ind)][:].split('\n'):
                    x = x.split()
                    if len(x) <2:
                        continue
                    pos.append(int(x[0]))
                    score.append(float(x[1]))
                del m[max(0,offset-self.L*(self.win+ind)) : offset+self.L*(self.win+1-ind)]

                # get the scores
                for i, base in enumerate(xrange(st, end)):
                    if base not in pos:
                        continue
                    self.score[num, i] = score[pos.index(base)]

                num += 1

        # moving average
        score = average(self.score, 0)
        self.mscore = zeros(score.shape,dtype=float)
        for x in range(self.ma, 2 * self.win + 1 - self.ma):
            self.mscore[x] = score[(x-self.ma): (x+self.ma+1)].mean()
        self.mscore = self.mscore[self.ma:self.ma*-1]



def MakeBinary(dir):
    '''MakeBinary(self, dir):
    '''
    for file in [open(x) for x in glob(dir+'/chr*')]:
        print file.name
        file.seek(-100, 2)
        data = zeros(int(file.read().split()[-2])+1,type=int16)
        file.seek(0)
        for line in file:
            pos, score = line.split()
            pos, score = int(pos), float(score) * 10000
            data[pos] = score
            if pos % 100000 == 0:
                print pos
        cPickle.dump(data, open(file.name + '.bin', 'w', 0), -1)
            
            
#read program options
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage = "%prog [options] bedfile list")
    
    parser.add_option("-d",  dest="dir", default = '/misc/iris/acct/weili/database/humanhg17_May2004/phastCons',
                      help="phastCons dir")
    
    parser.add_option("-r",  dest="range", default = 3000,type='int',
                      help="extend the whole region from the center to range, default 3000")

    parser.add_option("-w", dest='win', default=0,type='int',
                      help="bandwith for the moving average smoothing, default 0")
    
    parser.add_option("-m", dest='make', default=False,action="store_true",
                      help="make binary file")

    parser.add_option("--ymin",  dest="ymin", default = 0,type='float',
                      help="ymin in the plot")

    parser.add_option("--ymax",  dest="ymax", default = 0,type='float',
                      help="ymax in the plot")
    
    (options, args) = parser.parse_args()

    if len(sys.argv) <2:
        parser.print_help()
        sys.exit()
    
    if options.make:
        MakeBinary(sys.argv[1])
        sys.exit()
    
    for bed in args:
        out = bed.split('.')[0]
        phastCons(bed, options.dir, options.range, options.win, out, options.ymin, options.ymax)


'''
from rpy import r
self = phastCons('bedfile', '/misc/iris/acct/weili/database/humanhg17_May2004/phastCons', 3000, 250)
r.pdf('ConservationPlot.pdf')
r.plot(range(-3000, 3001), self.mscore, type = 'l', xlab = 'Distance from Center of Binding Sites', ylab = 'Conservation Score', lwd= 4, xlim = (-2000, 2000), ylim = (0.10, 0.25))
r.dev_off()
'''



