#!/usr/bin/env python
#wig 97-135374817


import os
import sys
import time

from optparse import OptionParser

def Info(infoStr):
    print "[%s] %s" %(time.strftime('%H:%M:%S'), infoStr)

# input a file with '\t' separate, return a 2-d list.
def TableFileInput(fileName, headDel = False, sep = "\t"):
    """TableFileInput(fileName, headDel = False)"""
    if type(fileName) == file:
        inf = fileName
    elif type(fileName) == str:
        if os.path.isfile(fileName):
            inf = open(fileName)
        else:
            print "no such file:<%s>" %fileName
            return 1
    l = []
    if headDel:
        inf.readline()
    for line in inf:
        line = line.rstrip("\n")
        l.append(line.split(sep))
    inf.close()
    return l

#file formated dataoutput (list/dict)
def TableFileOutput(l, filename, sep='\t'):
    """TableFileOutput(l,filename,RowNames='F')"""
    outf=file(filename,'w')
    for k in l:
        if type(k)==list:
            k=[str(m) for m in k]
            outf.writelines(sep.join(k)+'\n')
        else:
            outf.writelines(str(k)+'\n')
    outf.close()

def median(numbers):
    """Return the median of the list of numbers."""
    # Sort the list of numbers and take the middle element.
    n = len(numbers)
    copy = numbers[:] # So that "numbers" keeps its original order
    copy.sort()
    if n & 1: # There is an odd number of elements
        return copy[n/2]
    else:
        return (copy[n/2-1]+copy[n/2])/2

def prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first."""
    usage = "usage: %prog <-f prof_file -k key_file -m sort_method> [options]"
    description = "sort prof"
    # option processor
    optparser = OptionParser(version="%prog 0.0",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-f","--pfile",dest="pfile",type="string",help="input all the file you want to sort (also include the keyfile if you need). Multiple WIG files use ',' to split")
    optparser.add_option("-k","--kfile",dest="kfile",type="string",help="the ket file in sort file")
    optparser.add_option("-m",dest="method",type="string",help="sort method. choose from <median:median for each row>")
    optparser.add_option("-s",dest="suffix",type="string",help="the suffix for the output file")
    
    return optparser

def opt_validate(optparser):
    (options,args) = optparser.parse_args()
    options.pfile = options.pfile.split(',')
    if options.method not in ["median"]:
        Info("method can't find")
        sys.exit(1)
    pass
    return options

if __name__ == "__main__":
    opts=opt_validate(prepare_optparser())
    keylist = TableFileInput(opts.kfile, sep=',')
    width = len(keylist[0])
    index = iter(range(len(keylist)))
    t = [t.append(index.next()) for t in keylist]
    if opts.method == 'median':
        keylist.sort(key=lambda x:median([float(t) for t in x[:width]]))
    
    #read not key file
    for each in opts.pfile:
        pfile = TableFileInput(each, sep=',')
        out = []
        count = 0
        for i in keylist:
            out.append(pfile[i[-1]])
            count += 1
        if count != len(pfile):
            Info("Script Error.")
        TableFileOutput(out, "%s_%s" %(each, opts.suffix), sep=',')
