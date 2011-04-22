#!/usr/bin/env python
# weili@jimmy.harvard.edu



import random, math, os.path
import numpy.linalg as LA
from FileIO import *
import gaussian as G
from __builtin__ import round
from __builtin__ import sum as suml
from Wei import sample_wr, sd


class ProbeModel(Mybpmap):
    ''' ProbeModel inherits from FileIO.Mybpmap

        ##Methods:
        ProbeModel.CLEAN
        ProbeModel.DesignMatrix
        ProbeModel.FITMODEL
        ProbeModel.GetCels
        ProbeModel.MATCUTOFF
        ProbeModel.MATSCORE
        ProbeModel.MATscorePlot
        ProbeModel.MkUniqueIndex
        ProbeModel.NULLDIST
        ProbeModel.PRINTBED
        ProbeModel.READALL
        ProbeModel.REGIONCALL
        ProbeModel.Read
        ProbeModel.STANDARD

        ##Variables:
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
        self.outlier = 0                                        # consider outlier probes (True) or not (False)
        self.base = {'A':0,'C':1, 'G':2, 'T':3, 'a':0,'c':1, 'g':2, 't':3}
        self.affin = []                        # Affinity, the same shape as intensities
        self.beta = []                         # model parameter #cel * 81
        self.ndat = []                         # normalized(standard) data, the t(s)
        self.X = []                            # design matrix
        self.matscore = []                     # MAT score
        # MAT score for individual replicate, treat: tmean*sqrt(len) - pooled tmean_control*sqrt(len)
        # control: tmean_control*sqrt(len)
        self.matscore_m = []
        self.foldchange = []                   # trimmed mean alone, log fold-change
        self.MaxGap = ''                       # maximum gap between positive probes
        self.sigmahat = ''                     # Null distribution parameter
        self.mhat = ''                         # Null distribution parameter
        self.matcutoff = ''                    # interval cutoff, -10log10pvalue
        self.bed = []                          # interval regions
        self.nbed = {'P':1, 'N':1}             # number of bed regions, positive or negative
        self.ntrack = 1                        # track name
        self.fdrtable = []                     # different FDR at different MATscore cutoff
                                               # FDR(%)\t#NegPeak\t#PosPeak\tMATScore\t-10*log10_Pvalue
        # Affy  class
        self.seq = AffyFileParser.CGDACSequenceItem()
        self.hit = AffyFileParser.GDACSequenceHitItemType()
        self.bpmap = AffyFileParser.CBPMAPFileWriter()
    '''

    def __init__(self, bpmapname, GenomeGrp = 'Hs', outlier = False):
        Mybpmap.__init__(self, bpmapname, GenomeGrp)
        self.affin = []                        # Affinity, the same shape as intensities
        self.beta = []                         # model parameter #cel * 81
        self.ndat = []                         # normalized(standard) data, the t(s)
        self.X = []                            # design matrix
        self.matscore = []                     # MAT score tmean_treat*sqrt(len) - tmean_control*sqrt(leng)
        # MAT score for individual replicate, treat: tmean*sqrt(len) - pooled tmean_control*sqrt(len)
        # control: tmean_control*sqrt(len)
        self.matscore_m = []
        self.foldchange = []                   # trimmed mean alone, log fold-change
        self.MaxGap = ''                       # maximum gap between positive probes
        self.sigmahat = ''                     # Null distribution parameter
        self.mhat = ''                         # Null distribution parameter
        self.matcutoff = ''                    # interval cutoff, -10log10pvalue
        self.bed = []                          # interval regions
        self.nbed = {'P':1, 'N':1}             # number of bed regions, positive or negative
        self.ntrack = 1                        # track name
        self.fdrtable = []                     # different FDR at different MATscore cutoff
        self.outlier = outlier                 # consider outlier probes (True) or not (False)


    def DesignMatrix(self, PMProbe, copy):
        '''DesignMatrix(PMProbe, copy): make design matrix from PMProbe and copy_number
        PMProbe: probe sequence int8 array,
        copy:  copy number array
        Constructs a 81 pars X-matrix for the model
        3 * 25 "ACG" postions + 4 * 2 "ACGT" count and count square, 1 for copy number
        '''

        print >> sys.stderr, 'Making design matrix', time.asctime()
        X=zeros((PMProbe.shape[0], 81),dtype=float)
        # baseline
        X[:, 0] = sum(PMProbe == self.base['T'], 1).astype(float)
        j=1
        # 25 'ACG' postions
        for ibase in 'ACG':
            X[:, j:j+25] = PMProbe == self.base[ibase]
            j += 25
        # 'ACGT'  count square
        for ibase in 'ACGT':
            count = sum(PMProbe == self.base[ibase], 1).astype(float)
            X[:, j] = count**2
            j += 1
        # log copy number
        X[:, j] = log(copy) + 1
        return X


    def FITMODEL(self, n = 400000):
        '''FitModel( n = 400000): Fits the model with intensities and
        returns model parameter beta; unique n intensities and X; ExpIndex and UniqIndex, ReadIndex
        n[optional]: number of probes to fit, default 400000, set to 0 to fit all probes
        quick[optional]: 1: using evan"s quick least-sqare fitting;
        0: using scipy"s slow fitting with less memory
        Constructs a 81 pars X-matrix for the model from probe seqs and copy number
        1 for base line T-count, 3 * 25 "ACG" postions + 4  "ACGT"  count square, 1 for copy number
        '''

        self.MkUniqueIndex(nProbeRead = n)   # get full PMX, PMY, MatchScore
        self.Read(0, 'PMProbe')              # unique n PMProbe
        self.GetCels(self.PMX[self.ReadIndex], self.PMY[self.ReadIndex], self.outlier)   # unique n intensities
        # this probe is not an outlier in all replicates
        OutlierIndex = sum(self.intensities == 0, 1) == 0
        self.X = self.DesignMatrix(self.PMProbe, self.MatchScore[self.ReadIndex])  # ReadIndex X
        # release memory
        self.bpmap.Close()
        self.cels.Clean()
        for name in ['PMX', 'PMY', 'MatchScore', 'PMProbe']:
            self.__setattr__(name, [])
        print >> sys.stderr, 'Fitting model ... ', time.asctime()
        X = self.X[OutlierIndex]
        self.beta = dot(dot(LA.pinv(dot(transpose(X), X)), transpose(X)),
                        self.intensities[OutlierIndex])


    def READALL(self, N = 400000):
        ''' READALL(N = 400000): Get unique self.affin, self.intensities, and full self.Chr, self.Position,
        N [optional]: number of probes to get the affin at the same time, to save memory, default 400000
        '''

        # all the probes are used for model fitting
        if self.X.shape[0] == self.UniqIndex.shape[0] :
            print >> sys.stderr, 'Model fitting on all unique probes', time.asctime()
            self.affin = dot(self.X, self.beta)
            self.Read(1, 'Chr', 'Position')
            self.bpmap.Close()
        else:
            print >> sys.stderr, 'Too many probes on the array, \nrebuild the design matrix and get the affinity', \
                time.asctime()
            # re-read all the PMProbe data
            self.affin = zeros((self.UniqIndex.shape[0], self.beta.shape[1]),dtype=float)

            for i in range(0, self.UniqIndex.shape[0],  N):
                self.ReadIndex = self.UniqIndex[i: i+N].copy()
                print >> sys.stderr, self.ReadIndex
                self.Read(0, 'PMProbe', 'MatchScore')
                self.X = self.DesignMatrix(self.PMProbe, self.MatchScore)
                self.affin[i: i+N] = dot(self.X, self.beta)
            self.Read(1, 'Chr', 'Position', 'PMX', 'PMY')
            self.bpmap.Close()
            self.GetCels(self.PMX[self.UniqIndex], self.PMY[self.UniqIndex], self.outlier)
        # release memory
        self.cels.Clean()
        for name in ['PMProbe', 'MatchScore', 'X', 'PMX', 'PMY']:
            self.__setattr__(name, [])


    def STANDARD(self, binsize = 3000):
        '''STANDARD(binsize = 3000): returns model affinity centered and
        affinity-bin scaled data in self.ndat (the ts)
        binsize[optional]: number of affinity bins, 3000 probes for each bin by default
        '''

        self.ndat = zeros(self.affin.shape,dtype=float)
        # 3000 probes per affinity bin
        nUniq = self.UniqIndex.shape[0]
        n = nUniq / binsize

        for j in arange(self.ndat.shape[1]):
            print >> sys.stderr, 'Standardizing Sample:', os.path.basename(self.cels.names[j]), time.asctime()
            tmp = self.affin[:, j].argsort()
            for ilevel in arange(0.0, 1.0,  1.0/n):
                index = tmp[int(ilevel*nUniq) : int((ilevel+1.0/n) * nUniq)]
                if index.shape[0] < 5:
                    continue
                intens = self.intensities[index, j]
                affin = self.affin[index, j]
                binsds = intens[intens != 0].std()              # use non-outlier probes for bin SD
                self.ndat[index, j]  = (intens - affin) / binsds
        self.ndat[self.intensities == 0] = 0                   # outlier probe t-value 0
        # expand  ndat to full index
        self.ndat = self.ndat[self.ExpIndex]
        # write to disk
        self.cels.Write(str(self.outlier) + self.checksum, self.ndat)
        # release memory
        for name in ['intensities', 'affin']:
            self.__setattr__(name, [])



    def MATSCORE(self, group = '', BandWidth=300, short = 0, pair = '', var = 0, replicate = 0):
        '''MATSCORE(group = '', BandWidth=300, short = 0, pair = '', var = 0, replicate = 0):
        return self.matscore and/or self.matscore_m (The first nTreat columns are Treat data)
        BandWidth[optional]: BandWidth, default 300 bp, the total window will be 2*BandWidth
        short[optional]: minimum number of probes to get the trimmed mean, default is \
        half probes within the window. The matscore will be set to 0 if #probes < short.
        pair[optional]: subtract control from the treat at the t-value level
        var[optional]: control the input variance, divide sd(input) in MATScore, default no.
        replicate[optional]: set the replicate to 1 for matscore of each replicate self.matscore_m, default is 0 (no)
        '''

        # reading Position and Chr if necessary
        if len(self.Position) == 0:
            self.Read(1, 'Chr', 'Position')
            self.bpmap.Close()

        # reading .mat files to get the final self.ndat
        if self.cels.iscel.any():    # read the raw cel files
            ndat = zeros((self.cels.nintensity, self.Position.shape[0]),dtype=float)
            ndat[self.cels.iscel] = transpose(self.ndat)
            ndat[self.cels.iscel == False] = self.cels.ReadTvalue()
            self.ndat = ndat
        else:
            self.ndat = self.cels.ReadTvalue()

        # init
        print >> sys.stderr, "Making MAT score", time.asctime()
        print >> sys.stderr, "Control Input Variance : ", var
        self.trim = 0.1
        self.var = var
        start = time.clock()
        self.short = short
        if not short:
            self.short = int(round(BandWidth / 35.))
        self.BandWidth = BandWidth
        self.pair = pair
        # if not group, treat every cel as treatment
        if len(group) == 0:
            group = ones(self.ndat.shape[0],dtype=bool)
        self.group = group
        # subtract control from the treat at the t-value level
        if pair:
            self.ndat[group == True] = self.ndat[group == True] - self.ndat[group == False]

        nTreat = sum(group == True)
        nContr =  sum(group == False)
        Treat = ravel(transpose(self.ndat[group == True])).tolist()
        Control = ravel(transpose(self.ndat[group == False])).tolist()
        self.ndat.transpose()
        self.matscore = zeros(self.Position.shape,dtype=float)
        self.foldchange = zeros(self.Position.shape,dtype=float)
        matscoreTreat = zeros(self.Position.shape,dtype=float)
        if nContr > 0 and not self.pair:
            matscoreControl = zeros(self.Position.shape,dtype=float)
            ControlVar = ones(self.Position.shape,dtype=float)
        if replicate:
            self.matscore_m = zeros((self.Position.shape[0], nTreat + nContr),dtype=float)

        nhits = 0           # number of total hits
        for ichr, chr in enumerate(self.ChrInd):
            st, end = 0, 0                 # start, end of the sliding window [st, end) < win
            Position = self.Position[self.Chr == ichr]
            for i, pos in enumerate(Position):
                while((pos - Position[st]) > BandWidth):
                    st += 1
                while((Position[end] - pos) < BandWidth and end < Position.shape[0] - 1):
                    end += 1
                # minimum probe numbers
                if (end - st) < self.short:
                    continue

                # treats
                window = sorted(Treat[nTreat*(st+nhits): nTreat*(end+nhits)])
                bound = int(round(len(window) * self.trim))       # boundary for trimmed mean
                if  bound:
                    window = window[bound: -bound]
                mean_win = suml(window)/len(window)
                TreatSqrtlen = math.sqrt(len(window))
                matscoreTreat[i + nhits] = mean_win * TreatSqrtlen
                self.foldchange[i + nhits] = mean_win

                if replicate:
                    for irep in range(nTreat):
                        window = sorted(Treat[nTreat*(st+nhits) + irep : nTreat*(end+nhits) : nTreat])
                        bound = int(round(len(window) * self.trim))       # boundary for trimmed mean
                        if  bound:
                            window = window[bound: -bound]
                        mean_win = suml(window) / len(window)
                        self.matscore_m[i + nhits][irep] = mean_win * TreatSqrtlen  # use the pool treat sqrtlen

                # if controls and not pair
                if nContr > 0 and not self.pair:
                    window = sorted(Control[nContr*(st+nhits) : nContr*(end+nhits)])
                    bound = int(round(len(window) * self.trim))       # boundary for trimmed mean
                    if  bound:
                        window = window[bound: -bound]
                    mean_win = suml(window)/len(window)
                    sqrtlen = math.sqrt(len(window))
                    matscoreControl[i + nhits] = mean_win * sqrtlen
                    self.foldchange[i + nhits] -=  mean_win
                    # control the input variance,
                    if var:
                        ControlVar[i + nhits] = sd(window)

                    if replicate:
                        for irep in range(nContr):
                            window = sorted(Control[nContr*(st+nhits) + irep : nContr*(end+nhits) : nContr])
                            bound = int(round(len(window) * self.trim))
                            if  bound:
                                window = window[bound: -bound]
                            mean_win = suml(window)/len(window)
                            self.matscore_m[i + nhits][irep+ nTreat] = mean_win * TreatSqrtlen   # use the treat pool sqrtlen

                if i % 100000 == 0:
                    print >> sys.stderr, i, chr, time.clock() - start
            nhits += Position.shape[0]

        # re-scale matscore to get the same variance
        if nContr > 0 and not self.pair:
            ControSigmahat =  self.NULLDIST(matscoreControl)[1]
            ScaleFactor = self.NULLDIST(matscoreTreat)[1] /ControSigmahat
            ControlVar += 0.15 * ControSigmahat
            self.matscore = (matscoreTreat - matscoreControl * ScaleFactor)
            if var:
                self.matscore /= ControlVar
            self.mhat, self.sigmahat = self.NULLDIST(self.matscore)
            if replicate:
                matscoreControl.resize(matscoreControl.shape[0], 1)   # numpy returns None, not array copy
                self.matscore_m -= matscoreControl * ScaleFactor
                if var:
                    ControlVar.resize(ControlVar.shape[0], 1)         # numpy returns None, not array copy
                    self.matscore_m /= ControlVar
        else:
            self.matscore = matscoreTreat
            self.mhat, self.sigmahat = self.NULLDIST(self.matscore)



    def NULLDIST(self, dist):
        '''NullDist(dist): Get the Null distribution from left-half non-overlapping windows.
        return Null parameters, mhat, sigmahat
        '''
        nonz = dist != 0
        # non-overlapping matscore
        matscore = []
        for ichr, chr in enumerate(self.ChrInd):
            index = ((self.Chr == ichr)*nonz).astype(bool)
            Position = self.Position[index]
            if Position.shape[0] < 1:
                continue
            tt = dist[index]
            st = Position[0]
            for i, x in enumerate(Position):
                if x - st > self.BandWidth * 2:
                    matscore.append(tt[i])
                    st = x
        matscore = array(matscore)
        null = sort(matscore)[: matscore.shape[0]/2]
        mhat = null.max()
        sigmahat = concatenate((null, null * -1 + 2* mhat)).std()
        return mhat, sigmahat



    def MATCUTOFF(self,  FDR = 5,  Pvalue = '', Matscore = ''):
        '''MATCUTOFF(FDR = 5,  Pvalue = '', Matscore = ''): return the matscore cutoff
        FDR: we will find the maximum #(Pos-Neg)Peak in the fdrtable (<= FDR cutoff)
        and use the associated matscore as cutoff
        cutoff priority: FDR > Pvalue > Matscore
        '''
        if FDR != '':
            fdrtable = self.fdrtable[self.fdrtable[:,0] <= FDR]     # below this FDR cutoff
            self.matcutoff = fdrtable[fdrtable[:,3].argsort()[-1]][4]  # maximum real peaks
        elif Pvalue:
            zscore = G.gsl_cdf_ugaussian_Qinv(Pvalue)
            self.matcutoff =  zscore * self.sigmahat + self.mhat
        else:
            self.matcutoff = Matscore


    def FDRTABLE(self, MaxGap = 300, file = ''):
        '''FDRTABLE(self, MaxGap = 300, file = ''):
        return the FDR table
        FDR(%)	#NegPeak	#PosPeak	#(Pos-Neg)Peak	MATScore
        '''
        print >> sys.stderr, "Making FDR table", time.asctime()
        self.out = file
        if file:
            if self.ntrack == 1 :              # first array
                file = open(file+'.FDR.table.txt', 'w', 0)
            else:
                file = open(file+'.FDR.table.txt', 'a', 0)
        else:
            file = sys.stdout
        print >> file, '# ' + self.bpmapname
        print >> file, 'FDR(%)\t#NegPeak\t#PosPeak\t#(Pos-Neg)Peak\tMATScore\t-10*log10_Pvalue'
        self.fdrtable = []
        SMALL = 0.0001
        self.MaxGap = MaxGap
        NegMatscore = -1*self.matscore+2*self.mhat    # matscore for negative peak-call

        # start from pvalue 0.001
        cut = min(G.gsl_cdf_ugaussian_Qinv(0.001)* self.sigmahat + self.mhat,
                  self.matscore.max() - SMALL, NegMatscore.max() - SMALL)
        #cut = min(self.matcutoff,self.matscore.max() - SMALL,NegMatscore.max()-SMALL)
        #print "MATCUTOFF:",self.matcutoff
        #print "0001:",G.gsl_cdf_ugaussian_Qinv(0.001)* self.sigmahat + self.mhat
        # positive peaks
        pbed = array(self._BedCall(self.matscore, self.MaxGap, cut))[:,-1]
        # negative peaks
        nbed = array(self._BedCall(NegMatscore, self.MaxGap,  cut))[:,-1] + SMALL
        nbed = sort(nbed).tolist()
        nbed.reverse()

        for i, cut in enumerate(nbed):
            pvaluecut = -10 * log10(G.gsl_cdf_ugaussian_Q((cut - self.mhat)/self.sigmahat))
            if i < 30:
                pos = array(self._BedCall(self.matscore, self.MaxGap, cut)).shape[0]
                neg = array(self._BedCall(NegMatscore, self.MaxGap,  cut)).shape[0]
            else:
                pos = where(pbed >= cut)[0].shape[0]
                neg = i
            FDR = min(100.0,100.0 * neg / max(1, pos))
            self.fdrtable.append((FDR, neg, pos, pos - neg, cut, pvaluecut))
            print >> file, "%.2f\t%d\t%d\t%d\t%.2f\t%.2f" %  tuple(self.fdrtable[-1])
            if  FDR > 70 and pos > 100:
                break
        self.fdrtable = array(self.fdrtable)


    def FDRTABLEPlot(self, r = ''):
        FDR = self.fdrtable[:,0]
        Pos_Neg = self.fdrtable[:,3]
        r.plot(FDR, Pos_Neg, xlab = 'FDR%', ylab = 'Cummulative #Pos - Cummulative #Neg)', main = self.bpmapname+'\nFDR.table',
               pch = 20, cex = 1.5, ylim = (0, Pos_Neg.max()))


    def _BedCall(self, matscore, MaxGap, cut):
        '''_BedCall(self, matscore, MaxGap, cut):
        return bed list: chr st end maximum_scores
        '''
        self.bedhead = ['Chr', 'St', 'End', 'MATScore']
        bed = []
        # convert matscore to matrix
        matscore = matscore.copy()
        if len(matscore.shape) == 1:
            matscore.shape = matscore.shape[0], 1
        ind = sum((matscore >= cut), 1).astype(bool)
        chr = self.Chr[ind]
        pos = self.Position[ind]
        matscore = matscore[ind]

        # no peaks found
        if pos.shape[0] == 0 :
            return bed

        # coordinates are in the middle NOT at the beginning of probes
        ist = 0
        for i, ipos in enumerate(pos) :
            # last probe or next pos - current_pos > MaxGap or next_pos is in different chr
            if i == (pos.shape[0]-1) or \
            pos[i+1] - ipos >= MaxGap or \
            chr[i+1] != chr[i]:
                matscorebed = matscore[ist:i+1, :]
                bed1 = [chr[i], pos[ist]-12, ipos+12]
                # add the maximum matscore from each replicate within this region
                for irep in range(matscorebed.shape[1]):
                    bed1.append(matscorebed[:, irep].max())
                bed.append(bed1)
                ist = i+1
        return bed



    def _BedAnno(self, bed, extend = 0):
        '''_BedAnno(self, bed, extend = 0): Bed region annotation
        '''
        self.bedhead = ['Chr', 'St', 'End', 'Name', '-10log10_Pvalue', 'MATScore',
                        'Fold_Change', 'FDR(%)', 'Peak_Pos', 'Length']
        self.bedheadformat = ['%s', '%d', '%d', '%s', '%.2f', '%.2f',
                        '%.2f', '%.2f', '%d', '%d']
        if len(bed) == 0:
            return
        lastchr = ''
        for i, (chr, st, end) in enumerate(array(bed,dtype=int)[:,:3]):
            if chr != lastchr:                       # new chromosome
                ind = self.Chr == chr
                Position = self.Position[ind]
                matscore = self.matscore[ind]
                foldchange = self.foldchange[ind]
                lastchr = chr
                lastend = 1
            coordinates = where((Position >= st) * (Position <= end))

            # adding max-matscore-pvalue, peak_position, max_foldchange
            maxind = abs(matscore[coordinates]).argsort()[-1]
            max_matscore = matscore[coordinates][maxind]
            max_pvalue = max(1e-300, G.gsl_cdf_ugaussian_Q(abs(max_matscore - self.mhat)/self.sigmahat))
            max_pvalue = -10*log10(max_pvalue)
            max_fc = math.exp(foldchange[coordinates][maxind])
            FDR = self.fdrtable[self.fdrtable[:,4] <= max_matscore]
            if len(FDR):
                FDR = FDR[0][0]     # the highest one <  max_matscore
            else:                   # for negbed, negative matscore
                FDR = 100.0
            st = max(lastend + 1, (st - (extend + self.BandWidth)))         # never overlap with previous region
            end += (extend + self.BandWidth)
            fc_peak_position = Position[coordinates][maxind] - st
            # chr st end name max-matscore-pvalue max_fc fc_peak_position length
            bed[i] = [chr, st, end, '',  max_pvalue, max_matscore,  max_fc,
                           FDR, fc_peak_position,  end-st]
            lastend = end


    def REGIONCALL(self, cut = 0, MaxGap = 300, extend = 0):
        ''' RegionCall(cut = 0, MaxGap = 300, extend = 0): return ChIP-enriched regions
        The ChIP-region will be at least 2*BandWidth+1 long (if only one probe pass the cutoff)
        cut[optional]: default cutoff is self.matcutoff
        MaxGap[optional]: maximum gap between probes to merge probes into the same region, default 300 bp
        extend: bp to extend at both ends of the enriched region, default 0
        '''

        self.MaxGap = MaxGap
        if not cut:
            cut = self.matcutoff
        print >> sys.stderr, 'Region calling with cutoff', cut, time.asctime()
        self.bed = self._BedCall(self.matscore, MaxGap,  cut)
        self._BedAnno(self.bed, extend)
        # negative peaks
        self.negbed = self._BedCall(-1*self.matscore+2*self.mhat, MaxGap, cut)
        self._BedAnno(self.negbed, extend)



    def _PrintHeader(self, filehandle = sys.stdout):
        print >> filehandle, "#BandWidth=%d MaxGap=%d " % (self.BandWidth, self.MaxGap),
        if self.pair:
            print >> filehandle, 'Paired Test',
        if self.var:
            print >> filehandle, 'Control Input Variance ',
        if self.outlier:
            print >> filehandle, 'Remove outlier probes',
        print >> filehandle
        print >> filehandle, "#Null: mhat=%.2f sigmahat=%.2f" % (self.mhat, self.sigmahat)
        matscore = sort(self.matscore)
        print >> filehandle, "#MAT_Score: max=%.2f  min=%.2f  top_0.1%%=%.2f  bottom_0.1%%=%.2f" \
                % (matscore[-1], matscore[0], matscore[int(matscore.shape[0] * -0.001)],
                   matscore[int(matscore.shape[0] * 0.001)])
        print >> filehandle, "#Num_Pos_Peaks=%d  Num_Neg_Peaks=%d  FDR=%.2f%%" % \
                            (len(self.bed), len(self.negbed), 100.0*len(self.negbed)/max(1,len(self.bed)))
        print >> filehandle, "#bpmap:", self.bpmapname
        print >> filehandle, "#cel.Treat:", "   ".join([os.path.basename(x) for x in
                                                        self.cels.names[self.group]])
        print >> filehandle, "#cel.Control:", "   ".join([os.path.basename(x) for x in
                                                          self.cels.names[self.group == False]])


    def PRINTBED(self, file = '', blkname = 'Blk'):
        '''PRINTBED(file = '', blkname = 'Blk'):
        file[optional]: name of the output, default stdout
        print 1) .bed and .bed.xls files: ChIP-enriched regions
        2) .neg.xls file: Negative ChIP-enriched regions, i.e. swap the Treat and Control
        and call regions using exact the same cutoff  from 1)
        '''
        if file:
            if self.ntrack == 1 :              # clean the bed file
                open(file+'.bed', 'w', 0)
                print >> open(file+'.bed.xls', 'w', 0), \
                '%s\t'* len(self.bedhead) % tuple(self.bedhead)
                print >> open(file+'.neg.xls', 'w', 0),\
                '%s\t'* len(self.bedhead) % tuple(self.bedhead)
            bedfile = open(file+'.bed', 'a', 0)
            sumfile = open(file+'.bed.xls', 'a', 0)
            negfile = open(file+'.neg.xls', 'a', 0)
        else:
            negfile = bedfile = sumfile =  sys.stdout

        # no real peaks
        if len(self.bed) == 0:
            print >>bedfile, "track name=\"MATCh track%d\" description=\"No peaks\"" \
                                         % (self.ntrack)
            self._PrintHeader(sumfile)
            return

        # for browser position
        if self.ntrack == 1:
            chrind = [0] * (len(self.ChrInd))
            boundary = {}                               # chromosome boundaries
            for i, x in enumerate(self.bed):
                chr, st, end, name, score = x[0], x[1], x[2], x[3], x[4]
                chrind[chr] += 1
                # init
                if not boundary.has_key(chr):
                    boundary[chr] = [st, end]
                    continue
                if st < boundary[chr][0]:
                    boundary[chr][0] = st
                if end > boundary[chr][1]:
                    boundary[chr][1] = end
            chr = chrind.index(max(chrind))
            print >> bedfile, "browser position %s:%d-%d" % (self.ChrInd[chr], boundary[chr][0], boundary[chr][1])
            print >> bedfile, "track name=\"%s\"" % (file.split(os.path.sep)[-1])

        pvaluecut = G.gsl_cdf_ugaussian_Q((self.matcutoff - self.mhat)/self.sigmahat)
        FDR = 100.0*len(self.negbed)/len(self.bed)
        for filehandle in [sumfile, negfile]:
            print >>filehandle, "track name=\"%s_%d\" description=\"MATscore=%.2f Pvalue=%.2e FDR(%%)=%.2f\"" \
            % (file.split(os.path.sep)[-1], self.ntrack, self.matcutoff, pvaluecut, FDR)
        self._PrintHeader(sumfile)

        for ind in self.NumChr.argsort():
            for x in self.bed:
                if x[0] != ind:
                    continue
                for ifile in [bedfile, sumfile]:
                    print >> ifile, "%s\t%d\t%d\t%s_%d%s\t%.2f\t" % (self.ChrInd[x[0]], x[1], x[2], \
                                    blkname, self.nbed['P'], x[3], x[4]),
                print >> bedfile, ''
                print >> sumfile, '\t'.join(self.bedheadformat[5:]) % tuple(x[5:])
                self.nbed['P'] += 1
            # negative peaks
            for x in self.negbed:
                if x[0] != ind:
                    continue
                print >> negfile, "%s\t%d\t%d\t%s_%d%s\t%.2f" % (self.ChrInd[x[0]], x[1], x[2], \
                                    blkname, self.nbed['N'], x[3], x[4]),
                print >> negfile, '\t'.join(self.bedheadformat[5:]) % tuple(x[5:])
                self.nbed['N'] += 1
        for filehandle in [bedfile, sumfile, negfile]:
            print >> filehandle


    def PRINTTVALUE(self, file = ''):
        '''PRINTTVALUE(self, file = ''):
        '''
        print >> sys.stderr, 'printing t values ', file+'.'+ self.bpmapname +'.tsv'
        tfile = open(file +'.'+ self.bpmapname +'.tsv','w',0)
        print >> tfile, 'Chr\tPos\t','\t'.join([os.path.basename(x) for x in self.cels.names])
        for(x,y,z) in zip(self.ndat, self.Chr, self.Position):
            print >> tfile, str(self.ChrInd[y])+ '\t'+ str(z)+'\t'+ '\t'.join([str(i) for i in x])


    def HEATMAP(self, file='', r = ''):
        '''HEATMAP(self, file='', r = ''):
        '''
#        import numpy
        print >> sys.stderr, "HEATMAP, ", time.asctime()
        if not file:
            file = sys.stdout
        else:
            if self.ntrack == 1 :              # first track
                file = open(file+'.heatmap.txt', 'w', 0)
                print >> file, 'Chr\tSt\tEnd\tFinalMAT\t',
                cels = concatenate((self.cels.names[self.group], self.cels.names[self.group == False]))
                print >> file, '\t'.join([os.path.basename(x) for x in cels])
            else:
                file = open(file+'.heatmap.txt', 'a', 0)
            print >> file, '#', self.bpmapname
        if not self.matscore_m.max():
            self.MATSCORE(group = self.group, BandWidth=self.BandWidth, short = self.short,
                          pair = self.pair, var = self.var, replicate = True)
        matscore = self.matscore.copy()
        matscore.shape = matscore.shape[0], 1
        bed = array(self._BedCall(concatenate((matscore, self.matscore_m), 1),
                                  self.MaxGap, self.matcutoff))
        # no enriched regions
        if not len(bed):
            r.plot(9, main = 'No enriched regions for ' + self.bpmapname)
            return None
        # only one replicate
        if bed.shape[1] <= 5:
            r.plot(9, main = 'You need more than one replicate for this heatmap')
            return None

        negbed = array(self._BedCall(-1*concatenate((matscore, self.matscore_m),1),
                        self.MaxGap, self.matcutoff - 2*self.mhat))
        if len(negbed):
            negbed[:, 3:] *= -1
            bed = concatenate((bed, negbed), 0)
        for x in bed:
            print >> file, '%s\t%d\t%d\t' % (self.ChrInd[int(x[0])], x[1], x[2]),
            print >> file, '%.2f\t'*(len(x)-3) % tuple(x[3:])

        # no more than 6000 regions; save memory
        if bed.shape[0] > 6000:
            import random
            bed =  array(random.sample(bed, 6000))

        labCol = ['All']
        for x in range(self.group.sum()):
            labCol.append('T'+str(x+1))
        for x in range(sum(self.group == False)):
            labCol.append('C'+ str(x+1))

        pvaluecut = -10*log10(G.gsl_cdf_ugaussian_Q((self.matcutoff - self.mhat)/self.sigmahat))
        xlab = "#Neg #Pos #Call Matscore Pvalue\n%d %d %d %.2f %d" %(len(self.negbed),
                            len(self.bed),  bed.shape[0], self.matcutoff, pvaluecut)
        bed = bed[:, 3:]
        #r.heatmap(Rmatrix(bed), labCol=labCol, scale = 'none',
        #          main = self.bpmapname , xlab = xlab)
        # saturated at pos & neg cutoffs
        bed[bed < (-1*self.matcutoff + 2*self.mhat)] = (-1*self.matcutoff + 2*self.mhat)
        bed[bed > self.matcutoff] = self.matcutoff
        r.heatmap(r.as_matrix(array(bed)), labCol=labCol, scale = 'none',
                  main = self.bpmapname + ' [-1, 1]cutoff', xlab = xlab)
        return None


    def CLEAN(self):
        '''CLEAN(): clear associated data
        '''
        for name in ['affin', 'beta', 'ndat', 'X', 'intensities', 'matscore', 'pvalue', \
        'UniqIndex', 'ExpIndex', 'ReadIndex', 'Chr', 'MatchScore', 'PMProbe','PMX', 'PMY',\
        'Position', 'beta']:
            self.__setattr__(name, [])



    def MATscorePlot(self,  r = ''):
        ''' MATscorePlot(self, r = ''): MATScore histogram plot
        '''
        matscore = self.matscore[self.matscore != 0]
        nbreaks = matscore.shape[0]/400
        matscore_min = matscore.min() - 1
        r.par(mfrow=(2, 1))
        dd = r.hist(matscore, ylab = '', xlab = '', main = 'MAT score histogram ' + self.bpmapname,
                    breaks = nbreaks, border = 'red', col = 'red')
        dd = r.hist(matscore, ylab = '', xlab = '', main = 'MAT score histogram' + self.bpmapname,
                    breaks = nbreaks, border = 'red', col = 'red', ylim = (0, 25))
        r.par(mfrow=(1, 1))


    def ProbePlot(self,  r = ''):
        ''' ProbePlot(self, r = ''): probe intensity vs. copy number plot
        '''
        print >> sys.stderr, "probe intensity plot"
        # At probes
        GenomeGrp = self.GenomeGrp
        self.GenomeGrp = 'At'
        try:
            self.Read()
            self.GetCels()
            At_inten = self.intensities.copy()
        except:                             # no At probe
            At_inten = 0
        self.GenomeGrp = GenomeGrp
        self.CLEAN()

        self.MkUniqueIndex(nProbeRead = 400000)   # get full PMX, PMY, MatchScore
        self.bpmap.Close()
        self.GetCels(self.PMX[self.ReadIndex], self.PMY[self.ReadIndex])
        MatchScore = self.MatchScore[self.ReadIndex].round().astype(int)
        maxCopies = MatchScore.max()
        r.par(mfrow=(2, 4))
        for icel in range(self.intensities.shape[1]):
            names, boxplotdata = [], []
            if any(At_inten):
                names.append('At')
                boxplotdata.append(At_inten[:, icel])
            intensity = self.intensities[:, icel]
            for icpy in range(1, maxCopies):
                tmp = intensity[MatchScore == icpy][:10000]        # top 10k probes
                if tmp.shape[0]:
                    names.append(icpy)
                    boxplotdata.append(tmp)
            tmp = intensity[MatchScore > icpy]
            if tmp.shape[0]:
                names.append(icpy + 1)
                boxplotdata.append(tmp)
            r.boxplot(boxplotdata, names = names, ylab = 'Log Probe Intensity',
                   xlab = 'Copy Number', pch = 19, col = 'red',
                   main = os.path.basename(self.cels.names[icel]))
            print >> sys.stderr, os.path.basename(self.cels.names[icel])
        r.par(mfrow=(1, 1))


    def ProbeTvaluePlot(self, r = ''):
        '''probe t-value (fitting without copy number) vs copy number plot
        '''
        r.par(mfrow=(2, 4))
        print >> sys.stderr, "probe t-value plot"
        # re-read everything
        if not len(self.ReadIndex) or not (len(self.ReadIndex) == len(self.intensities)
                                        == len(self.MatchScore)):
            self.MkUniqueIndex(nProbeRead = 400000)   # get full PMX, PMY, MatchScore
            self.GetCels(self.PMX[self.ReadIndex], self.PMY[self.ReadIndex])
            self.MatchScore = self.MatchScore[self.ReadIndex]
        self.Read(0, 'PMProbe')
        intensities = self.intensities
        MatchScore = self.MatchScore.round().astype(int)
        maxCopies = MatchScore.max()
        X = self.DesignMatrix(self.PMProbe, ones(self.PMProbe.shape[0]))
        OneCpy_X = X[MatchScore == 1]
        beta = dot(dot(LA.pinv(dot(transpose(OneCpy_X), OneCpy_X)), transpose(OneCpy_X)),
                   intensities[MatchScore == 1])

        # At probes
        try:
            GenomeGrp = self.GenomeGrp
            self.GenomeGrp = 'At'
            self.Read()
            self.GetCels()
            intensities = concatenate((intensities, self.intensities))
            X = concatenate((X, self.DesignMatrix(self.PMProbe, ones(self.PMProbe.shape[0]))))
            MatchScore = concatenate((MatchScore, zeros(self.PMProbe.shape[0])))
            self.GenomeGrp = GenomeGrp
        except:
            pass
        affin = dot(X, beta)
        nUniq = affin.shape[0]
        binsize = 3000
        n = nUniq / binsize
        for icel in range(intensities.shape[1]):
            tmp = affin[:, icel].argsort()
            ndat = zeros(nUniq,dtype=float)
            for ilevel in arange(0.0, 1.0, 1.0/n):
                index = tmp[int(ilevel*nUniq) : int((ilevel+1.0/n) * nUniq)]
                if index.shape[0] < 5:
                    continue
                intens = intensities[index, icel]
                aff = affin[index, icel]
                binsds = intens[intens != 0].std()              # use non-outlier probes for bin SD
                ndat[index]  = (intens - aff) / binsds
            names, boxplotdata = [], []
            # At probes
            tmp = ndat[MatchScore == 0]
            if len(tmp):
                names.append('At')
                boxplotdata.append(tmp)
            for icpy in range(1, maxCopies):
                tmp = ndat[MatchScore == icpy][:10000]        # top 10k probes
                if tmp.shape[0]:
                    names.append(icpy)
                    boxplotdata.append(tmp)
            tmp = ndat[MatchScore > icpy]
            if tmp.shape[0]:
                names.append(icpy + 1)
                boxplotdata.append(tmp)
            self.box = boxplotdata
            r.boxplot(boxplotdata, names = names, ylab = 'Probe t-value',
                   xlab = 'Copy Number', pch = 19, col = 'red',
                   main = os.path.basename(self.cels.names[icel]),
                   ylim = (max(-10, ndat.min()), min(10, ndat.max())))
            print >> sys.stderr, os.path.basename(self.cels.names[icel])
        r.par(mfrow=(1, 1))



    def __ReadBed(self, bedfile):
        ''' testing
        '''

        self.bed = []
        for x in open(bedfile):
            x = x.split()
            if x[0][:3] != 'chr':
                continue
            if x[0] not in self.ChrInd:
                continue
            ind = self.ChrInd.index(x[0])
            Position = self.Position[self.Chr == ind]
            st, end = int(x[1]), int(x[2])
            coordinates = where((Position >= st) * (Position <= end) )
            if coordinates[0].shape[0] < 9:
                continue
            self.bed.append([ind, st,end])
