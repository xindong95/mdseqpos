#!/usr/bin/env python
# weili@jimmy.harvard.edu

import os.path, sys, copy, re
from numpy import *

'''
change space to center-to-center distance
'''

class Bed(object):
    ''' Bed(self, filename, space = 0, percent = 0): Bed file object
    space[default 0]: two regions are considered the same if the gap between two regions is
    < gap. gap of -200 means that they have to have at least 200 bp overlap to be consided as the same.
    percent[default 0]: two regions are considered the same if they overlap more than
    this percent(based on the Bed object itself). For now, you can set up to 50% overlap.
    Return:
    self.data[chr] = [st, end, score,.. , name_index]
    self.name = names

    ### code ####
    for x in self.data[chr]:
        st, end, score, name_id =  x[0], x[1], x[2:-1], int(x[-1])
        print st, end, score, self.name[name_id]
    '''

    def __init__(self, filename, space = 0, percent = 0.5):
        self.filename = os.path.basename(filename)
        self.data = {}                              # data[chr] = [st, end, score,.. , name_index]
        self.name = []                              # name str array
        self.space = space                          # distance for overlapping regions
        self.percent = percent
        rawdata = []
        for x in open(filename):
            x = x.lower()
            if re.search(r'chr\S', x[:4], re.I):
                x = x.split()
                rawdata.append(x)
        if len(rawdata) :
            rawdata = char.array(rawdata)

            numlist = ones(rawdata.shape[1],dtype=bool)      # boolean list, 1 for numeric 0 otherwise
            for i, x in enumerate(rawdata[0]):
                try:
                    x = float(x)
                except:
                    numlist[i] = 0

            for chr in set(rawdata[:,0]):
                self.data[chr] = []
            numdata = transpose(transpose(rawdata)[numlist])
            ns = numdata.shape
            numdata = array([float(x) for x in ravel(numdata)]) # numpy returns None not array copy upon resize
            numdata.resize(ns)
            for i, x in enumerate(numdata.tolist()):
                self.data[rawdata[i, 0]].append(x + [i])                    # add line index as the last element
            self.name = char.array(transpose(transpose(rawdata)[numlist == 0][1:]), itemsize = 100)
            if not len(self.name):                              # chr as the name
                self.name = rawdata[:,0:1]
            self._sort()


    def _sort(self):
        ''' Sort bed regions
        '''
        for chr in self.data.keys():
            self.data[chr] = array(self.data[chr])
            ind = self.data[chr][:,0].argsort()
            self.data[chr] = self.data[chr][ind]


    def __repr__(self):
        out = []
        id = 0                  # number of regions
        chrs = sorted(self.data.keys())
        for chr in chrs:
            for x in self.data[chr]:
                st, end, score, name_id =  x[0], x[1], x[2:-1], int(x[-1])
                id += 1
                out.append("%s\t%d\t%d\t" % (chr, st, end) + '\t'.join(self.name[name_id]) +
                           '\t' + '\t'.join(["%.2f" % (x) for x in score]))
        if id > 20:
            out = char.array(out)
            return "\n".join(out[:10]) + '\n...\n' + '\n'.join(out[-10:]) + '\n#entries %d'%(id)
        else:
            return "\n".join(out) + '\n#entries %d'%(id)


    def __str__(self):
        out = []
        id = 0
        chrs = sorted(self.data.keys())
        for chr in chrs:
            for x in self.data[chr]:
                st, end, score, name_id =  x[0], x[1], x[2:-1], int(x[-1])
                id += 1
                out.append("%s\t%d\t%d\t" % (chr, st, end) + '\t'.join(self.name[name_id]) +
                           '\t' + '\t'.join(["%.2f" % (x) for x in score]))
        return "\n".join(out)


    def _newcoord(self, st, end):
        percentgap = int((end - st) * self.percent )
        middle = int(0.5 * (st + end))
        st = min(middle, st - self.space + percentgap)
        end = max(middle, end + self.space - percentgap)
        return st, end


    def __sub__(self, bed2):
        ''' remove all peaks in self, which are falling within bed2
        '''
        bed = copy.deepcopy(self)
        for chr in bed.data.keys():
            if bed2.data.has_key(chr):
                for i, (st, end) in enumerate(bed.data[chr][:, :2]):
                    st, end = self._newcoord(st, end)
                    # overlap
                    if ((bed2.data[chr][:, 0] - end)  * (bed2.data[chr][:,1] - st)).min() <= 0:
                        bed.data[chr][i] = 0
                bed.data[chr] =  bed.data[chr][(bed.data[chr] != 0)[:,0]]
        for chr in bed.data.keys():
            if bed.data[chr].shape[0] == 0:
                bed.data.pop(chr)
        return bed


    def __add__(self, bed2):
        ''' union of bed and bed2
        '''
        bed = copy.deepcopy(self)
        bed2 = copy.deepcopy(bed2)
        name_len = bed.name.shape[1]
        score_len = bed.data.values()[0].shape[1]
        name2_len = bed2.name.shape[1]
        score2_len = bed2.data.values()[0].shape[1]

        # expand bed.name to include names from bed2
        addon = char.array([' '] * bed.name.shape[0] * name2_len, shape = (bed.name.shape[0],
                                                name2_len))
        bed.name = concatenate((bed.name, addon),1)
        for chr in bed.data.keys():
            data =  bed.data[chr]
            # expand bed.data to include scores from bed2
            bed.data[chr] = concatenate((data[:,:-1], zeros((data.shape[0],score2_len -1)),
                                         data[:, -1:]), 1)
            if bed2.data.has_key(chr):
                data2 = bed2.data[chr]
                for i, (st, end) in enumerate(data[:, :2]):
                    st, end = self._newcoord(st, end)
                    offset = (data2[:, 0] - end)  * (data2[:,1] - st)
                    #  overlap
                    if offset.min() <= 0:
                        overlap_id = offset.argsort()[0]          # closest overlap
                        bed.data[chr][i][score_len -1: -1] = data2[overlap_id][:-1]
                        bed.name[int(data[i][-1])][-1* name2_len :] = \
                                bed2.name[int(data2[overlap_id][-1])]
                        bed2.data[chr][overlap_id] = 0          # replace overlap entry in bed2 as 0

        num_name = bed.name.shape[0]
        bed.name = bed.name.tolist()
        # bed2-alone region
        addon = char.array([' '] * bed2.name.shape[0] * name_len, shape = (bed2.name.shape[0],
                                                name_len))
        bed2.name = concatenate((addon, bed2.name),1)
        for chr in bed2.data.keys():
            data = bed2.data[chr]
            data = data[data[:,0] != 0]  # remove overlap entries in bed2
            if len(data) == 0:
                continue
            # expand bed2.data
            bed2.data[chr] = concatenate((data[:,:2],zeros((data.shape[0], score_len-1)), data[:,2:]),1)
            if not bed.data.has_key(chr):
                bed.data[chr] = []
            else:
                bed.data[chr] = bed.data[chr].tolist()
            for x in bed2.data[chr]:
                bed.name.append(bed2.name[int(x[-1])])    #x[-1]: name_index
                x[-1] = num_name
                bed.data[chr].append(x)
                num_name += 1
            bed.data[chr] = array(bed.data[chr])
        bed.name = char.array(bed.name)
        bed._sort()
        return bed


    def __and__(self, bed2):
        ''' intersection of bed and bed2
        '''
        bed = copy.deepcopy(self)
        name_len = bed.name.shape[1]
        score_len = bed.data.values()[0].shape[1]
        name2_len = bed2.name.shape[1]
        score2_len = bed2.data.values()[0].shape[1]

        # expand bed.name to include names from bed2
        addon = char.array([' '] * bed.name.shape[0] * name2_len, shape = (bed.name.shape[0],
                                                name2_len))
        bed.name = concatenate((bed.name, addon),1)
        for chr in bed.data.keys():
            data =  bed.data[chr]
            # expand bed.data to include scores from bed2
            bed.data[chr] = concatenate((data[:,:-1], zeros((data.shape[0],score2_len -1)), \
                                         data[:, -1:]), 1)
            if bed2.data.has_key(chr):
                data2 = bed2.data[chr]
                for i, (st, end) in enumerate(data[:, :2]):
                    st, end = self._newcoord(st, end)
                    offset = (data2[:, 0] - end)  * (data2[:,1] - st)
                    #  overlap
                    if offset.min() <= 0:
                        overlap_id = offset.argsort()[0]
                        bed.data[chr][i][score_len -1: -1] = data2[overlap_id][:-1]
                        bed.name[int(data[i][-1])][-1* name2_len :] = \
                                bed2.name[int(data2[overlap_id][-1])]
                    else:
                        bed.data[chr][i] = 0
                bed.data[chr] = bed.data[chr][(bed.data[chr] != 0)[:,0]]
                if not bed.data[chr].shape[0]:
                    bed.data.pop(chr)
            else:
                bed.data.pop(chr)
        return bed


    def __or__(self, bed2):
        ''' only return bed itself
        '''
        bed = copy.deepcopy(self)
        name_len = bed.name.shape[1]
        score_len = bed.data.values()[0].shape[1]
        name2_len = bed2.name.shape[1]
        score2_len = bed2.data.values()[0].shape[1]

        # expand bed.name to include names from bed2
        addon = char.array([' '] * bed.name.shape[0] * name2_len, shape = (bed.name.shape[0],
                                                name2_len))
        bed.name = concatenate((bed.name, addon),1)
        for chr in bed.data.keys():
            data =  bed.data[chr]
            # expand bed.data to include scores from bed2
            bed.data[chr] = concatenate((data[:,:-1], zeros((data.shape[0],score2_len -1)), \
                                         data[:, -1:]), 1)
            if bed2.data.has_key(chr):
                data2 = bed2.data[chr]
                for i, (st, end) in enumerate(data[:, :2]):
                    st, end = self._newcoord(st, end)
                    offset = (data2[:, 0] - end)  * (data2[:,1] - st)
                    #  overlap
                    if offset.min() <= 0:
                        overlap_id = offset.argsort()[0]
                        bed.data[chr][i][score_len -1: -1] = data2[overlap_id][:-1]
                        bed.name[int(data[i][-1])][-1* name2_len :] = \
                                bed2.name[int(data2[overlap_id][-1])]
        return bed


    def __getitem__(self, k):
        # init
        bed = copy.deepcopy(self)
        bed.data = {}
        out = []
        id = 0
        chrs = self.data.keys()
        chrs.sort()
        if type(k) == NumArray:
            k = k.tolist()

        if type(k) == int:
            for chr in chrs:
                for i, x in enumerate(self.data[chr]):
                    if id == k:
                        bed.data[chr] = self.data[chr][i:i+1]
                        return bed
                    id += 1

        elif type(k) == list:
            for chr in chrs:
                for x in self.data[chr]:
                    if k.count(id):                         # id in k
                        if not bed.data.has_key(chr):
                            bed.data[chr] = []
                        for i in range(k.count(id)):
                            bed.data[chr].append(x)
                    id += 1
                if bed.data.has_key(chr):
                    bed.data[chr] = array(bed.data[chr])
            return bed

        elif type(k) == slice :      # return bed object
            if not k.start:
                k = slice(0, k.stop, k.step)
            if not k.stop:
                k = slice(k.start, 1e99, k.step)
            for chr in chrs:
                for x in self.data[chr]:
                    if id >= k.start and id <  k.stop:
                        if not bed.data.has_key(chr):
                            bed.data[chr] = []
                        bed.data[chr].append(x)
                    id += 1
                if bed.data.has_key(chr):
                    bed.data[chr] = array(bed.data[chr])
            return bed

        elif type(k) == tuple:      # return array object
            k  = list(k)
            if not k[0].start:
                k[0] = slice(0, k[0].stop, k[0].step)
            if not k[0].stop:
                k[0] = slice(k[0].start, 1e99, k[0].step)
            for chr in chrs:
                for item in self.data[chr]:
                    if id >= k[0].start and id < k[0].stop:
                        if k[1] == 0:               # 0 for chr
                            out.append(chr)
                        elif k[1] > 0:              # positive number for other numeric columns
                            out.append(item[k[1]-1])
                        else:                       # negative for name
                            out.append(self.name[int(item[-1])][ int(-1 * k[1] -1) ])
                    id += 1
            if max(k[1], 0) == 0:
                out = char.array(out)
            else:
                out = array(out)
            return out

        else:
            return k


    def __len__(self):
        id = 0
        chrs = self.data.keys()
        chrs.sort()
        for chr in chrs:
            id += self.data[chr].shape[0]
        return id


    def RmRepeat(self, *arg):
        bed = copy.deepcopy(self)
        pattern = []
        for x in arg:
            pattern.append(re.compile(x))
        for chr in bed.data.keys():
            for i, ind in enumerate(bed.data[chr][:,-1]):
                name = bed.name[int(ind),0]
                for x in pattern:
                    if x.search(name):
                        bed.data[chr][i] = 0
                        print >> sys.stderr,  name
                        break
            bed.data[chr] =  bed.data[chr][(bed.data[chr] != 0)[:,0]]
            if bed.data[chr].shape[0] == 0:
                bed.data.pop(chr)
        return bed


    def Merge(self, gap):
        '''Merge(self, gap): merge bed regions, regions separated by < gap will be merged into one
        '''
        bed = copy.deepcopy(self)
        for chr in bed.data.keys():
            data =  bed.data[chr].copy()
            st = -1e9999
            for ipos , pos in enumerate(data[:, 0]):
                if pos - st < gap:
                    # replace the current st to the previous st
                    data[ipos, 0] = data[ipos-1, 0]
                    data[ipos-1] = 0
                st = pos
            bed.data[chr] = data[(data != 0)[:,0]]
        return bed

    def Encode2Bed(self):
        '''Encode2Bed(self): Transfer the encode bed format (chr pos score) to real bed format
        '''
        bed = copy.deepcopy(self)
        chrs = sorted(self.data.keys())
        for chr in chrs:
            d = bed.data[chr]
            bed.data[chr] = concatenate((d[:,0:1], d[:,0:1],  d[:,1:]), 1)
        return bed


    def Change(self, ColId, newitem):
        '''Change(self, ColId, newitem): index is the column in the self.data[chr] matrix,
        we will replace that column with the newitem
        '''
        bed = copy.deepcopy(self)
        id = 0
        chrs = self.data.keys()
        chrs.sort()

        for chr in chrs:
            for i in range(len(bed.data[chr])):
                bed.data[chr][i, ColId] = newitem[id]
                id += 1
        return bed

    def ChangePeak(self):
        '''Change the coordinates to the peak position
        '''
        peak = self[:,7]+self[:,1]
        return self.Change(0, peak).Change(1, peak)


    def ChangeMiddlePeak(self):
        '''Change the coordinates to be the middle of the peak position and the physical center
        '''
        peak = (0.5* (self[:,1]+ self[:,7] + 0.5 * (self[:,1]+self[:,2]))).astype(Int32)
        return self.Change(0, peak).Change(1, peak)

    def Top(self, n = 100, ColId = 4):
        '''Top(self, n = 100, ColId = 4):
        return top n (default 100) sites
        '''
        bed = copy.deepcopy(self)
        return bed[bed[:,ColId].argsort()[-n:]]
