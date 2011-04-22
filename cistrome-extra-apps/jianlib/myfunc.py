import time, os, sys

# print current time and information on screen
def Info(infoStr):
    sys.stderr.write("[%s] %s\n" %(time.strftime('%H:%M:%S'), infoStr))


# input a file with '\t' separate, return a 2-d list.
def TableFileInput(fileName, headDel = False, sep = "\t"):
    """TableFileInput(fileName, headDel = False)"""
    if type(fileName) == file:
        inf = fileName
    elif type(fileName) == str:
        if os.path.isfile(fileName):
            inf = open(fileName)
        else:
            Info("no such file:<%s>" %fileName)
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

def Orderfile(keyfile, pfile, kmethod, sep=','): #pfile is a list
    keylist = TableFileInput(keyfile, sep=',')
    #print keylist[4726:4729]
    width = len(keylist[0])
    #print width
    index = iter(range(len(keylist)))
    map(lambda x:x.append(index.next()), keylist)
    #print keylist[4726:4729]
    #t = [t.append(index.next()) for t in keylist]
    if kmethod == "median":
        keylist.sort(key=lambda x:median([float(t) for t in x[:width]]))
    elif kmethod == "maximum":
        keylist.sort(key=lambda x:max([float(t) for t in x[:width]]))
    elif kmethod == "mean":
        keylist.sort(key=lambda x:sum([float(t) for t in x[:width]]))
    
    #read not key file
    for each in pfile:
        pfile = TableFileInput(each, sep=',')
        out = []
        count = 0
        for i in keylist:
            out.append(pfile[i[-1]])
            count += 1
        if count != len(pfile):
            Info("Sth Error.")
            sys.stderr.write("%d != %d\n" %(count, len(pfile)))
            sys.stderr.write("pfile=%s\n" %each)
            
        TableFileOutput(out, each, sep=',')

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


