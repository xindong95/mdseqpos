#!/usr/bin/env python
# weili@jimmy.harvard.edu


import sys, time, os.path, AffyFileParser, re, ConfigParser, zipfile, os, md5
import cPickle as pickle
from numpy import *
from glob import glob



class Zip(object):
    def __init__(self, filename):
        files = glob(filename+'.*.bar')
        zip = zipfile.ZipFile(filename + '.bar.zip', 'w', zipfile.ZIP_DEFLATED)
        dirname = os.path.dirname(files[0])
        if dirname:
            os.chdir(dirname)
        for f in files:
            zip.write(os.path.basename(f))
        zip.close()


class Tag(object):
    '''Tag(filename): TAG (Tiling Array Group) file parser
    self.cel = {1:[a.cel,b.cel], 2:[c.cel,d.cel]}
    self.bpmap = {1:a.bpmap, 2:b.bpmap}
    '''

    def __init__(self, filename):
        self.cel = {}
        self.bpmap = {}
        parser = ConfigParser.ConfigParser()
        self.parser = parser

        #############################
        # obligatory options
        ##############################
        try:
            if not parser.read(filename):
                raise Exception, ('Cannot find tag file '+ filename)
            self.filename = re.compile('.tag', re.I).sub('', os.path.basename(filename))
            self.dir = os.path.dirname(filename)
            if self.dir :
                self.dir += os.path.sep

            # output
            self.out = self.dir + self.filename
            if parser.get('output', 'log'):
                sys.stdout = sys.stderr = open(self.out +'.log', 'a', 0)
            print >> sys.stderr, '[ ' + time.asctime() + ' ]'

            # data
            for x in ['bpmapfolder', 'celfolder']:
                if parser.get('data', x)[-1] != os.path.sep:
                    parser.set('data',x, parser.get('data',x)+os.path.sep)
            self.bpmapfolder = parser.get('data', 'bpmapfolder')
            self.celfolder = parser.get('data', 'celfolder')
            if not os.access(self.celfolder, os.W_OK):
                raise Exception, ('No write access to the cel folder %s. Please either change the directory access permission or hard link(or copy) your cel files to a writable directory' %(self.celfolder))

            self.genomegrp = parser.get('data', 'genomegrp')

            # cel
            if not parser.options('cel') or not parser.options('bpmap'):
                raise Exception, 'No cel or bpmap files'

            # group; 1=treatment, 0=control, "#"=skip
            if parser.get('data', 'group'):
                self.group = array([int(x) for x in parser.get('data', 'group') if x.isdigit()],dtype=bool)
                group = array([x.isdigit() for x in parser.get('data', 'group')],dtype=bool)
            else:
                self.group = ones(len(parser.items('cel')[0][1].split()),dtype=bool)
                group = ones(len(self.group),dtype=bool)

            # check cels and bpmaps
            for i in sorted(parser.options('bpmap')):
                self.bpmap[i] = parser.get('bpmap', i)
                tmp = AffyFileParser.CBPMAPFileWriter()
                tmp.SetFileName(self.bpmapfolder + self.bpmap[i])
                if not tmp.ReadHeader():
                    raise Exception, ('Not a valid bpmap file ' + tmp.GetFileName())
                self.cel[i] = char.array(parser.get('cel', i).split())
                if len(self.cel[i]) != len(group):
                    raise Exception, ('The number of cel files in ' + self.bpmap[i] + ' is not the same as #elements in the [data] -> Group parameter. Please separate all your cel files into Treatment (1) or Control (0) groups (e.g. 111000 for 3 ChIP and 3 Input, 10 for 1 ChIP and 1 Input, etc). To skip a cel file, use #, e.g.11##00 for 2 ChIP, 2 files to skip and 2 Input.')
                
                self.cel[i] = self.cel[i][group]
                for cel in self.cel[i]:
                    tmp = AffyFileParser.CCELFileData()
                    tmp.SetFileName(self.celfolder + cel)
                    if not tmp.ReadHeader():
                        raise Exception, ('Not a valid cel file ' + tmp.GetFileName())

                # printing
                print >> sys.stderr, '\n', self.bpmap[i]
                print >> sys.stderr, 'Treat: ', ' '.join(self.cel[i][self.group])
                print >> sys.stderr, 'Control: ', ' '.join(self.cel[i][self.group == 0])

            # intensity analysis
            self.bandwidth = parser.getint('intensity analysis', 'bandwidth')
            if self.bandwidth <= 0:
                raise Exception, 'BandWidth error, should be positive integer'
            self.maxgap = parser.getint('intensity analysis', 'maxgap')
            if self.maxgap <= 0:
                raise Exception, 'MaxGap error, should be positive integer'
            self.minprobe = parser.getint('intensity analysis', 'minprobe')
            if self.minprobe < 0 :
                raise Exception, 'MinProbe error, should positive integer'

            # interval analysis
            if parser.get('interval analysis', 'pvalue'):
                self.pvalue = parser.getfloat('interval analysis', 'pvalue')
                self.fdr = self.matscore = ''
                if self.pvalue <= 0:
                    raise Exception, 'Error: Pvalue cutoff should > 0'
            elif parser.get('interval analysis', 'matscore'):
                self.matscore = parser.getfloat('interval analysis', 'matscore')
                self.pvalue = self.fdr = ''
            else:
                self.fdr = 0.01
                self.matscore = self.pvalue = ''

        except:
            print >> sys.stderr, 'Tag file format error: ',  sys.exc_info()[0], sys.exc_info()[1]
            sys.exit()
        #############################
        # misc options (not necessary)
        ##############################

        # repeat library
        try:
            self.replib = parser.get('data', 'replib')
        except:
            self.replib = ''
        if self.replib:
            try:
                open(self.replib, 'rb')               # file test
            except:
                print >> sys.stderr, 'Cannot open repeat library ', self.replib
                sys.exit()

        # could be matscore or foldchange or both
        try:
            self.barfile = parser.get('misc', 'barfile').split()
        except:
            self.barfile = ['matscore']

        # retrieve genomic sequences
        try:
            self.genome = parser.get('misc', 'genome')
        except:
            self.genome = ''

        # map bed regions to genome annotation
        try:
            self.annotation = parser.get('misc', 'annotation').split()
            self.annotationfolder = parser.get('misc', 'annotationfolder')
            if self.annotationfolder[-1] != os.path.sep:
                self.annotationfolder += os.path.sep
        except:
            self.annotation = self.annotationfolder = ''

        # get the paired information
        try:
            self.pair = parser.getboolean('data', 'pair')
            if self.pair and where(self.group == False)[0].shape[0] != where(self.group)[0].shape[0]:
                raise Exception, '#Treatments and #Controls are not the same, no pairing for this project.'
        except:
            self.pair = False

        # extend the bed region
        try:
            self.extend = parser.getint('interval analysis', 'extend')
        except:
            self.extend = 0

        # compress the bar files
        try:
            self.zip = parser.getboolean('misc', 'zip')
        except:
            self.zip = False

        # phastCons score
        try:
            self.phastCons_dir = parser.get('misc', 'phastCons_dir')
        except:
            self.phastCons_dir = ''
        try:
            self.phastCons_range = parser.getint('misc', 'phastCons_range')
            self.phastCons_win = parser.getint('misc', 'phastCons_win')
        except:
            self.phastCons_range = 3000
            self.phastCons_win = 250

        # control the variance from the input, default NO
        try:
            self.var = parser.getint('intensity analysis', 'var')
        except:
            self.var = 0

        # generate a .tsv file for raw t-value of every probe, default False
        try:
            self.tvalue = parser.getboolean('intensity analysis', 'tvalue')
        except:
            self.tvalue = False

        # FDR cutoff
        try:
            self.fdr = parser.getfloat('interval analysis', 'FDR')
            self.matscore = self.pvalue = ''
            if self.fdr < 0 or self.fdr > 100:
                raise Exception, 'Error: FDR cutoff should be within [0, 100]'
        except:
            pass

        # diagnoses
        try:
            self.diag = parser.getboolean('interval analysis', 'Diag')
        except:
            self.diag = False

        # remover outlier probes, default False
        try:
            self.outlier = parser.getboolean('intensity analysis', 'Outlier')
        except:
            self.outlier = False

        # output wiggle file instead of bar file
        try:
            self.wig = parser.getboolean('output', 'Wig')
        except:
            self.wig = False



class Mybpmap(object):
    '''Mybpmap(filename, GenomeGrp = "Hs"): read GenomeGroup data in bpmapfile
    GenomeGrp[optional]: name of the genomic group in this array, 'Hs' for Affy
    Human Tiling 1.0 and 2.0, ''(blank) for Affy Chr21/22 and Encode tiling array
    '''

    def __init__(self, bpmapname, GenomeGrp = 'Hs'):
        '''SetFileName(bpmapname, GenomeGrp = 'Hs'): Set Filename and GenomeGroup
        '''
        self.Chr = []                                           # chromosome information
        self.Position = []                                      # positions
        self.MatchScore = []                                    # copy number
        self.PMProbe = []                                       # probe seqs
        self.PMX = []                                           # PMX
        self.PMY = []                                           # PMY
        self.MMX = []                                           # MMX
        self.MMY = []                                           # MMY
        self.ChrInd = []                                        # index of chromosome information
        self.NumChr = []                                        # index of numerical chromosome
        self.pars = ['Chr','MatchScore', 'Position',  'PMProbe', 'PMX', 'PMY']   # elements in the array
        self.intensities = []                                   # cel intensities
        self.UniqIndex = []                                     # index for unique probes
        self.ExpIndex = []                                      # X[UniqIndex][ExpIndex] =  X
        self.ReadIndex = []                                     # index for reading probes
        self.base = {'A':0,'C':1, 'G':2, 'T':3, 'a':0,'c':1, 'g':2, 't':3}
        self.checksum = md5.new(open(bpmapname).read()).hexdigest()     # md5 checksum for the bpmap file
        # Affy  class
        self.seq = AffyFileParser.CGDACSequenceItem()
        self.hit = AffyFileParser.GDACSequenceHitItemType()
        self.bpmap = AffyFileParser.CBPMAPFileWriter()

        self.GenomeGrp = GenomeGrp                              # name of the genomic group in this array
        self.bpmap.SetFileName(bpmapname)
        self.bpmapname = os.path.basename(bpmapname)
        if not self.bpmap.Exists():
            print >> sys.stderr, 'Cannot find ' , bpmapname
            sys.exit()
        print >> sys.stderr, "Reading", self.bpmapname,  time.asctime()



    def Read(self,  all = 1, *pars):
        '''Read(all = 1, *pars):  Read array elements
        all[optional]: if True read all elements, if False, read elements based on self.ReadIndex
        *pars[optional]: can be "Chr","MatchScore"(copy_number/1e6), "Position",
        "PMProbe", "MMProbe", "PMX", "PMY, MMX, MMY"
        So Read() will read all elements of all the probes on the array
        Read(0, "PMProbe") will read PMProbe only based on self.ReadIndex
        Read(1, "PMProbe") will read all PMProbe
        result will be self.PMProbe, self.MatcshScore, etc.
        self.Chr will be index array from ChrInd
        self.PMProbe will be int8 two-dimentional array with {"A":0,"C":1, "G":2, "T":3}
        '''

        # init
        if not pars:
            pars = self.pars
        else:
            pars = list(pars)
        print >> sys.stderr,  '  '.join(pars), time.asctime()
        if all:
            print >> sys.stderr, 'All probes'
        else:
            print >> sys.stderr, 'Partial probes', self.ReadIndex.shape[0]
        if self.bpmap.GetNumberSequences() == 0:
            if not self.bpmap.Read():
                print >> sys.stderr, 'Error: ', self.bpmap.GetError()
                sys.exit()
        nhits = 0                                   # number of total hits
        for name in pars:
            self.__setattr__(name, [])
        if 'MMProbe' in pars:
            self.PMProbe = []
        self.ChrInd = []

        # read
        for iseq in range(self.bpmap.GetNumberSequences()):
            self.bpmap.GetSequenceItem(iseq, self.seq)
            if self.seq.GroupName() != self.GenomeGrp:
                continue
            print >> sys.stderr, 'reading ', self.seq.GetName()
            self.ChrInd.append(self.seq.GetName())
            ChrInd = len(self.ChrInd) - 1
            # index for this seq
            if not all:
                index = self.ReadIndex[((self.ReadIndex < (nhits + self.seq.GetNumberHits()))\
                    * (self.ReadIndex >= nhits)).astype(bool)]
                index = index - nhits
            else:
                index = arange(self.seq.GetNumberHits())
            for ihit in index:
                self.seq.GetHitItem(ihit, self.hit, 1)
                for name in pars:
                    if name == 'Chr':
                        self.Chr.append(ChrInd)
                    elif name in ['PMProbe', 'MMProbe']:
                        self.PMProbe.append([self.base[ibase] for
                            ibase in self.hit.PMProbe])
                    else:
                        self.__getattribute__(name).append(self.hit.__getattribute__(name))
            nhits += self.seq.GetNumberHits()
        if len(self.ChrInd) == 0:
            raise Exception, ('No data for this GenomeGrp ', self.GenomeGrp)
        # make array
        for name in ['Position', 'PMX', 'PMY', 'MMX', 'MMY']:
            if pars.count(name) > 0:
                self.__setattr__(name, array(self.__getattribute__(name),dtype=int))
        if pars.count('PMProbe') > 0:
            self.PMProbe = array(self.PMProbe,dtype=int8)
        if pars.count('MMProbe') > 0:
            self.PMProbe = array(self.PMProbe,dtype=int8)
            self.PMProbe[:, 12] = 3 - self.PMProbe[:, 12]
        if pars.count('Chr') > 0:
            self.Chr = array(self.Chr,dtype=int8)
        # MatchScore * 1e6 = copy number
        if pars.count('MatchScore') > 0:
            self.MatchScore = array(self.MatchScore, dtype=float) * 1e6

        # sort chromosomes correctly; chrM = 100001, chrX = 100002, chrY = 100003
        self.NumChr = []
        MaxNumChr = 100000
        for x in self.ChrInd:
            if x[3:].isdigit():
                self.NumChr.append(int(x[3:]))
            else:
                MaxNumChr += 1
                self.NumChr.append(MaxNumChr)
        self.NumChr = array(self.NumChr)


    def MkUniqueIndex(self, nProbeRead = 400000, ReUniq = 1):
        '''MkUniqueIndex(nProbeRead = 400000): return ReadIndex for unique probes according to
        PMX, PMY, MatchScore. return self.UniqIndex, self.ExpIndex, self.ReadIndex
        nProbeRead[optional]: number of probes to read, default 400000, set it to 0 to read all probes
        ReUniq[optional]: re-calculate the UniqIndex, ExpIndex, default True
        '''

        if ReUniq:
            if not len(self.PMX):
                self.Read(1, 'PMX', 'PMY', 'MatchScore')
            print >> sys.stderr, 'Making Unique Index ',  time.asctime()
            print >> sys.stderr, 'Maximum copy number: ', int(round(self.MatchScore.max())),\
            ' duplicate probe measurements: ', (self.MatchScore > 1).astype(int).sum()
            index = arange(0, self.PMX.shape[0])
            dpos = index[self.MatchScore > 1]            # list of index of duplicated probes
            # PMX, PMY of duplicated probes
            PMX = self.PMX[self.MatchScore > 1]
            PMY = self.PMY[self.MatchScore > 1]
            dprobes = array([(x<<32) + PMY[i] for i, x in enumerate(PMX)])
            dic = {}        # probe -> index
            for i, x in enumerate(dprobes):
                if not dic.has_key(x):
                    dic[x] = []
                dic[x].append(dpos[i])
            # the first occurance of dup probes
            unique = array([x[0] for x in dic.values()])
            # put the first occurance of dup probes back
            self.UniqIndex = sort(concatenate((index[self.MatchScore <= 1], unique)))

            #  X[UniqIndex][ExpIndex] =  X
            self.ExpIndex = zeros(self.PMX.shape,dtype=int)
            undic = {}
            for x in dic.values():
                undic[x[0]] = x
            for i, x in enumerate(self.UniqIndex):
                if undic.has_key(x):
                    self.ExpIndex[undic[x]] = i
                else:
                    self.ExpIndex[x] = i

        if nProbeRead:
            ratio =  max(1, 1.0 * self.UniqIndex.shape[0] / nProbeRead)
        else:
            ratio = 1
        nUniq = self.UniqIndex.shape[0]
        self.ReadIndex = array([self.UniqIndex[i] for i in arange(0, nUniq, ratio,dtype=int) if i < nUniq])


    def GetCels(self, PMX = [], PMY = [], outlier = False):
        '''GetCels( PMX = [], PMY = [], outlier = False): return the cel intensities, log_PM, from self.cels
        if not (PMX and PMY), will use self.PMX, self.PMY instead
        outlier[optional]: assign 0 to outlier probes, default NO.
        '''

        print >> sys.stderr, 'Getting cel intensities: ',  time.asctime()
        # no cel reading
        if self.cels.cels[0].GetNumCells() == 0:
            self.cels.Read()
        if (len(PMX) == 0 or len(PMY) == 0):
            if not len(self.PMX):
                self.Read(1, 'PMX', 'PMY')
            PMX = self.PMX
            PMY = self.PMY
        self.intensities = [self.cels.GetIntensity(X, PMY[i], outlier) for i, X in enumerate(PMX)]
        print >> sys.stderr, 'finished self.cels.GetIntensity'
        self.intensities =  log(array(self.intensities)).astype(float)


    def GetMMCels(self):
        '''GetMMCels(self): return the cel intensities, log_MM, from self.cels
        '''

        print >> sys.stderr, 'Getting cel intensities: ',  time.asctime()
        # no self.cels.cels
        if not len(self.cels.cels):
            for name in self.cels.names:
                cel = AffyFileParser.CCELFileData()
                cel.SetFileName(name)
                if not cel.Exists():
                    raise Exception, ('Not a valid cel file ', name)
                self.cels.cels.append(cel)
        # no cel reading
        if self.cels.cels[0].GetNumCells() == 0:
            self.cels.Read()
        if not len(self.MMX):
            self.Read(1, 'MMX', 'MMY')
        PMX = self.MMX
        PMY = self.MMY
        self.intensities = [self.cels.GetIntensity(X, PMY[i]) for i, X in enumerate(PMX)]
        self.intensities =  log(array(self.intensities)).astype(float)



class Mycels(object):
    '''Mycels([celfiles_list], head = 0): return cel object list as self.cels
    head is sef.outlier + md5.hexdigest on the bpmap files,
    if the bpmap file is not the original one used to
    generate the .mat files, we will discard the .mat files and re-read the cel files
    also returned the t-values associated with each cel file as self.ndat
    '''

    def __init__(self, names, head = 0):
        self.iscel = []    # 0 for t-values in .mat file, 1 for raw cel file
        self.cels = []
        self.nintensity = len(names)
        for name in names:
            tname = self.__tname(name)
            # read .mat file directly
            if os.path.isfile(tname) and pickle.load(open(tname))[0] == head:
                self.iscel.append(False)
            else:
                # read cel file
                self.iscel.append(True)
                cel = AffyFileParser.CCELFileData()
                cel.SetFileName(name)
                if not cel.Exists():
                    raise Exception, ('Not a valid cel file ', name)
                self.cels.append(cel)
        self.iscel = array(self.iscel,dtype=bool)
        self.names = char.array(names)


    def Read(self):
        '''Read(): Read cel files
        '''
        for  cel in self.cels:
            print >> sys.stderr, 'reading ', os.path.basename(cel.GetFileName()), time.asctime()
            if not cel.Read():
                print >> sys.stderr, cel.GetError(), os.path.basename(cel.GetFileName())
                sys.exit()


    def ReadTvalue(self):
        '''ReadTvalue(): Read t-values, .mat files
        '''
        if self.iscel.all():  # no t-values files to read
            return None
        ndat = []
        for x in where(self.iscel == False)[0]:
            tname = self.__tname(self.names[x])
            print >> sys.stderr, 'reading ', os.path.basename(tname), time.asctime()
            ndat.append(pickle.load(open(tname))[1])
        return array(ndat, dtype=float)


    def GetIntensity(self, x, y, outlier = False):
        ''' getintensity( x, y,  outlier=False)
        x, y: position of the probe on the array
        outlier[optional]: assign value 1 to outlier probes (logPM = 0)
        '''
        intensity = []
        for cel in self.cels:
            if outlier and cel.IsOutlier(x, y):
                intensity.append(1)
            else:
                intensity.append(cel.GetIntensity(x, y))
            if intensity[-1] <= 0:                    # cel file check
                raise Exception, ('Abnormal value %s at X %s Y %s, please check your cel file %s' % (
                    intensity[-1], x, y, cel.GetFileName()))
        return intensity


    def Clean(self):
        '''Clean(): Close all cel files to release memory
        '''
        for cel in self.cels:
            cel.Clear()


    def Write(self, head, ndat):
        '''Write(head, ndat): Write binary file of the ndat to .mat files
        '''
        for i, cel in enumerate(self.cels):
            name = self.__tname(cel.GetFileName())
            pickle.dump((head, ndat[:,i].copy()), open(name, 'wb',0), -1)


    def __tname(self, celname):
        pattern = re.compile('.CEL$', re.IGNORECASE)   # cel file suffix pattern
        if pattern.search(celname):
            return pattern.sub('.mat', celname)
        else:
            return celname + '.mat'



class Mybars(object):
    '''Save and Read bar files
    '''

    def __init__(self):
        self.bar = AffyFileParser.CBARFileWriter()
        self.data = AffyFileParser.BarSequenceResultData()
        self.seq = AffyFileParser.CGDACSequenceResultItem()


    def Save(self, bp, filename = '', content = ['matscore']):
        '''Save(bp, filename, , content = ["matscore"]): Save bpmap file
        bp: ProbeModel class instance
        filename: bar file name head, the final bar file will be filename.bpmapname_pvalue(intensity).bar
        content: could be ["matscore", "foldchange", "tvalue", "PMProbe"]
        '''

        print >> sys.stderr, "Saving bar files", time.asctime()
        filename = filename+'.'+os.path.basename(bp.bpmapname)
        for name in content:
            if name == 'tvalue':
                self.SaveTvalue(bp, filename)
                continue
            if name == 'PMProbe':
                self.SavePMProbe(bp, filename)
                continue
            if name == 'MMProbe':
                self.SaveMMProbe(bp, filename)
                continue

            # for 'matscore', 'foldchange'
            value = bp.__getattribute__(name)
            if name == 'matscore':
                self.bar.SetFileName(filename+'_matscore.bar')
                self.bar.CreateNewFile()
                self.bar.AddAlgorithmParameter('file_type', 'matscore')
            elif name == 'foldchange':
                self.bar.SetFileName(filename+'_foldchange.bar')
                self.bar.CreateNewFile()
                self.bar.AddAlgorithmParameter('file_type', 'FoldChange')

            self.bar.AddAlgorithmParameter('scale', 'NaturalLog')
            self.bar.AddAlgorithmParameter('genomic_map', os.path.basename(bp.bpmapname))
            self.bar.AddAlgorithmParameter('BandWidth', str(bp.BandWidth))
            self.bar.AddAlgorithmParameter('mhat', str(bp.mhat))
            self.bar.AddAlgorithmParameter('sigmahat', str(bp.sigmahat))
            self.bar.AddAlgorithmParameter('var', str(bp.var))
            self.bar.AddAlgorithmParameter('MaxGap', str(bp.MaxGap))
            self.bar.AddAlgorithmParameter('GenomeGrp', str(bp.GenomeGrp))
            self.bar.AddAlgorithmParameter('MinProbe', str(bp.short))
            # cel.treat
            for icel, cel in enumerate(bp.cels.names[bp.group]):
                self.bar.AddAlgorithmParameter('TreatCel_'+str(icel), os.path.basename(cel))
            # cel.control
            for icel, cel in enumerate(bp.cels.names[bp.group == False]):
                self.bar.AddAlgorithmParameter('ControlCel_'+str(icel), os.path.basename(cel))
            self.bar.AddColumn(AffyFileParser.BAR_DATA_INTEGER)
            self.bar.AddColumn(AffyFileParser.BAR_DATA_FLOAT)
            self.bar.SetNumberSequences(len(bp.ChrInd) + 1)

            nhits = 0                                           # total number of hits
            for i in range(self.bar.GetNumberSequences()-1):
                p = self.bar.GetResultsPtr(i)
                id = where(bp.Chr == i)[0]
                p.SetNumberDataPoints(id.shape[0])
                p.SetName(bp.ChrInd[i])
                p.SetGroupName('')
                p.SetVersion('')
                for j in range(p.GetNumberDataPoints()):
                    self.data.iValue = bp.Position[j + nhits]
                    p.SetDataPoint(j, 0, self.data)
                    self.data.fValue =  value[j + nhits]
                    p.SetDataPoint(j, 1, self.data)
                nhits += p.GetNumberDataPoints()

            self.bar.Save()
            self.bar.Close()

    def SaveWiggle(self, bp, filename = '', content = ['matscore']):
        '''Save(bp, filename, , content = ["matscore"]): Save bpmap file
        bp: ProbeModel class instance
        filename: bar file name head, the final bar file will be filename.bpmapname_pvalue(intensity).bar
        content: could be ["matscore", "foldchange", "tvalue", "PMProbe"]
        '''

        print >> sys.stderr, "Saving Wiggle files", time.asctime()
        filename = filename+'.'+os.path.basename(bp.bpmapname)+".MATscore.wig"
        fhd = open(filename,"w")
        fhd.write("track type=wiggle_0 name=\"MAT score for %s\" description=\"MAT score for %s\"\n" % (filename,filename))
        for name in content:
            if name == 'tvalue':
                self.SaveTvalueWiggle(bp, filename)
                continue

            # for 'matscore', 'foldchange'
            value = bp.__getattribute__(name)

            # cel.treat
            #for icel, cel in enumerate(bp.cels.names[bp.group]):
            #    self.bar.AddAlgorithmParameter('TreatCel_'+str(icel), os.path.basename(cel))
            # cel.control
            #for icel, cel in enumerate(bp.cels.names[bp.group == False]):
            #    self.bar.AddAlgorithmParameter('ControlCel_'+str(icel), os.path.basename(cel))
            #self.bar.AddColumn(AffyFileParser.BAR_DATA_INTEGER)
            #self.bar.AddColumn(AffyFileParser.BAR_DATA_FLOAT)
            self.bar.CreateNewFile()
            self.bar.SetNumberSequences(len(bp.ChrInd) + 1)

            nhits = 0                                           # total number of hits
            for i in range(self.bar.GetNumberSequences()-1):
                p = self.bar.GetResultsPtr(i)
                id = where(bp.Chr == i)[0]
                p.SetNumberDataPoints(id.shape[0])
                p.SetName(bp.ChrInd[i])
                fhd.write("variableStep chrom=%s span=25\n" % bp.ChrInd[i])
                for j in range(p.GetNumberDataPoints()):
                    fhd.write("%d\t%f\n" % (bp.Position[j + nhits], value[j + nhits]))
                nhits += p.GetNumberDataPoints()

        fhd.close()

    def SaveTvalue(self, bp, filename = '', content = 'tvalue'):
        '''SaveTvalue(self, bp, filename = ''):
        '''
        for icel, cel  in enumerate(bp.cels.names):
            celname = os.path.basename(cel)
            value = bp.ndat[:,icel]
            self.bar.SetFileName(filename+'.'+celname+'.tvalue.bar')
            self.bar.CreateNewFile()
            self.bar.AddAlgorithmParameter('file_type', 'tvalue')
            self.bar.AddColumn(AffyFileParser.BAR_DATA_INTEGER)
            self.bar.AddColumn(AffyFileParser.BAR_DATA_FLOAT)
            self.bar.SetNumberSequences(len(bp.ChrInd) + 1)
            nhits = 0                                           # total number of hits
            for i in range(self.bar.GetNumberSequences()-1):
                p = self.bar.GetResultsPtr(i)
                id = where(bp.Chr == i)[0]
                p.SetNumberDataPoints(id.shape[0])
                p.SetName(bp.ChrInd[i])
                p.SetGroupName('')
                p.SetVersion('')
                for j in range(p.GetNumberDataPoints()):
                    self.data.iValue = bp.Position[j + nhits]
                    p.SetDataPoint(j, 0, self.data)
                    self.data.fValue =  value[j + nhits]
                    p.SetDataPoint(j, 1, self.data)
                nhits += p.GetNumberDataPoints()
            self.bar.Save()
            self.bar.Close()

    def SaveTvalueWiggle(self, bp, filename = '', content = 'tvalue'):
        '''SaveTvalue(self, bp, filename = ''):
        '''
        for icel, cel  in enumerate(bp.cels.names):
            celname = os.path.basename(cel)
            fhd = open(filename+'.'+celname+'.tvalue.wig',"w")
            fhd.write("track type=wiggle_0 name=\"MAT score for %s\" description=\"MAT score for %s\"\n" % (filename,filename))
            value = bp.ndat[:,icel]
            self.bar.CreateNewFile()
            self.bar.SetNumberSequences(len(bp.ChrInd) + 1)
            nhits = 0                                           # total number of hits
            for i in range(self.bar.GetNumberSequences()-1):
                p = self.bar.GetResultsPtr(i)
                id = where(bp.Chr == i)[0]
                p.SetNumberDataPoints(id.shape[0])
                p.SetName(bp.ChrInd[i])

                fhd.write("variableStep chrom=%s span=25\n" % bp.ChrInd[i])
                for j in range(p.GetNumberDataPoints()):
                    fhd.write("%d\t%f\n" % (bp.Position[j + nhits], value[j + nhits]))

                nhits += p.GetNumberDataPoints()
            fhd.close()


    def SavePMProbe(self, bp, filename = ''):
        '''SavePMProbe(self, bp, filename = ''):
        '''
        bp.GetCels()
        bp.cels.Clean()
        # scale to the same mean
        # celmean = sum(bp.intensities, 0) / bp.intensities.shape[0]
        # bp.intensities = bp.intensities / celmean * 10.0
        for icel, cel  in enumerate(bp.cels.names):
            celname = os.path.basename(cel)
            value = bp.intensities[:,icel]
            self.bar.SetFileName(filename+'.'+celname+'.logPM.bar')
            self.bar.CreateNewFile()
            self.bar.AddAlgorithmParameter('file_type', 'logPM')
            self.bar.AddColumn(AffyFileParser.BAR_DATA_INTEGER)
            self.bar.AddColumn(AffyFileParser.BAR_DATA_FLOAT)
            self.bar.SetNumberSequences(len(bp.ChrInd) + 1)
            nhits = 0                                           # total number of hits
            for i in range(self.bar.GetNumberSequences()-1):
                p = self.bar.GetResultsPtr(i)
                id = where(bp.Chr == i)[0]
                p.SetNumberDataPoints(id.shape[0])
                p.SetName(bp.ChrInd[i])
                p.SetGroupName('')
                p.SetVersion('')
                for j in range(p.GetNumberDataPoints()):
                    self.data.iValue = bp.Position[j + nhits]
                    p.SetDataPoint(j, 0, self.data)
                    self.data.fValue =  value[j + nhits]
                    p.SetDataPoint(j, 1, self.data)
                nhits += p.GetNumberDataPoints()
            self.bar.Save()
            self.bar.Close()


    def SaveMMProbe(self, bp, filename = ''):
        '''SavePMProbe(self, bp, filename = ''):
        '''
        bp.GetMMCels()
        bp.cels.Clean()
        # scale to the same mean
        # celmean = sum(bp.intensities, 0) / bp.intensities.shape[0]
        # bp.intensities = bp.intensities / celmean * 10.0
        for icel, cel  in enumerate(bp.cels.names):
            celname = os.path.basename(cel)
            value = bp.intensities[:,icel]
            self.bar.SetFileName(filename+'.'+celname+'.logMM.bar')
            self.bar.CreateNewFile()
            self.bar.AddAlgorithmParameter('file_type', 'logMM')
            self.bar.AddColumn(AffyFileParser.BAR_DATA_INTEGER)
            self.bar.AddColumn(AffyFileParser.BAR_DATA_FLOAT)
            self.bar.SetNumberSequences(len(bp.ChrInd) + 1)
            nhits = 0                                           # total number of hits
            for i in range(self.bar.GetNumberSequences()-1):
                p = self.bar.GetResultsPtr(i)
                id = where(bp.Chr == i)[0]
                p.SetNumberDataPoints(id.shape[0])
                p.SetName(bp.ChrInd[i])
                p.SetGroupName('')
                p.SetVersion('')
                for j in range(p.GetNumberDataPoints()):
                    self.data.iValue = bp.Position[j + nhits]
                    p.SetDataPoint(j, 0, self.data)
                    self.data.fValue =  value[j + nhits]
                    p.SetDataPoint(j, 1, self.data)
                nhits += p.GetNumberDataPoints()
            self.bar.Save()
            self.bar.Close()


    def Read(self, bp, filename):
        '''Read(bp, filename): Read the intensity (trimmed meas score) into bp object.
        bp: ProbeModel class instance
        filename: bar file name head, the final bar file will be filename.bpmapname_intensity.bar
        '''

        filename = filename+'.'+os.path.basename(bp.bpmapname)+ '_signal.bar'
        print >> sys.stderr, "Reading bar files", filename, time.asctime()
        self.bar.SetFileName(filename)
        if not self.bar.Read():
            print >> sys.stderr, self.bar.GetError()
            sys.exit()

        # clean data in bp object
        for x in ['Chr', 'Position', 'matscore', 'spvalue',]:
            bp.__setattr__(x, [])
        bp.ChrInd = []

        # file header
        if self.bar.GetParameter(0).Value != 'signal':
            print >> sys.stderr, 'Not a signal bar file, please check'
            sys.exit()
        bp.bpmapname = self.bar.GetParameter(2).Value
        bp.BandWidth = int(self.bar.GetParameter(3).Value)
        bp.mhat = float(self.bar.GetParameter(5).Value)
        bp.sigmahat = float(self.bar.GetParameter(6).Value)

        for i in range(self.bar.GetNumberSequences() - 1):
            self.bar.GetResults(i, self.seq)
            bp.ChrInd.append(self.seq.GetName())
            for j in range(self.seq.GetNumberDataPoints()):
                self.seq.GetData(j, 0, self.data)
                bp.Position.append(self.data.iValue)
                self.seq.GetData(j, 1, self.data)
                bp.matscore.append(self.data.fValue)
                bp.Chr.append(i)

        # make array
        for x in ['Position', 'matscore', 'spvalue', 'Chr']:
            bp.__setattr__(x, array(bp.__getattribute__(x)))
