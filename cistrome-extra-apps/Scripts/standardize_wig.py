#!/usr/bin/env python

import string
import sys

def findSpan(file):
    """ Takes in a specific wig file and returns the span of the file in integer form. """
    f=open(file,'r')
    spanTmp=''
    text = f.read()
    spanPos = text.find('span')
    i=5
    while text[spanPos+i]!='\n':
        spanTmp=spanTmp+text[spanPos+i]
        i+=1
    span=int(spanTmp)
    return span

def chromDict(datalines,span):
    """ Takes in a list of the datalines of a wig file and the span of the wig file.
    RETURNS:
    datadict- a dictionary containing keys are values which corresepond to the posistion and relative signals given in the wig file.  The values given in the dict are for one chromosome
    stopLine- the stopLine is the heading line the describes the next chromosome in the wig file
    datalines- a new list of the datalines of the wig file without the lines that were used to make the dataDict
    lastNumb- gives the value of the last basepair that is described in dataDict."""
    dataDict={}
    i=0
    stopLine=''
    lastNumb = None
    for l in datalines:
        if l.startswith('track') or l.startswith('browse'):
            i+=1
            continue
        line=l.replace('\n','')
        dataArray=string.split(line)
        try:
            location = int(dataArray[0])
        except ValueError:
            stopLine = l
            datalines = datalines[i+1:len(datalines)]
            break
        value = dataArray[1]
        for x in range(location,location+span):
            dataDict[x]=value
        lastNumb = location
        i+=1
    return [dataDict, stopLine, datalines, lastNumb]
        
def findChromNumb(file):
    """ Takes in a wig file and returns the number of chromosomes the file describes in integer form. """
    f = open(file,'r')
    text = f.read()
    chromNumb = text.count('span')
    return chromNumb

def changeSpan(line,interval):
    """ Takes in a declaration line such as, 'variableStep chrom=chr1 span=10 and a new interval and returns a declaration line with the span number changed to the new interval."""
    spanPos = line.find('span')
    tmpLine= line[0:spanPos+5]
    newLine=tmpLine+str(interval)
    return newLine

def standardizeWig(file,interval,newFileName):
    """ Takes in a wig file and creates a new wig file with the new desired interval.  The file name of the new file is defined by newFileName. """
    newFile = open(newFileName,"w")
    f=open(file,'r')
    datalines=f.readlines()
    for line in datalines:
        if line.startswith('track') or line.startswith('browse'):
            newFile.write(line)
        else: break
    span = findSpan(file)
    chromNumb = findChromNumb(file)
    for x in range(chromNumb+1):
        print datalines
        [dataDict,stopLine,datalines,lastNum]=chromDict(datalines,span)

        i=1
        if lastNum != None:
            while i<=lastNum+span:
                if i in dataDict:
                    newFile.write(str(i)+'\t'+str(dataDict[i])+'\n')
                i+=interval
        if stopLine!='':
            newStopLine = changeSpan(stopLine,interval)
            newFile.write(newStopLine+'\n')
    return     
        
if __name__== '__main__':
    if len(sys.argv) == 4:
        standardizeWig(sys.argv[1], int(sys.argv[2]), sys.argv[3])
    else:
        print """USAGE:
        python standardize_wig.py [src_wig_filename] [interval] [output_filename]
Example:
        python standardize_wig.py test.wig 128 test_128.wig
        """


