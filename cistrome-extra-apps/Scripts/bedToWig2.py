#!/usr/bin/env python

import string
import sys
import operator

def linInterpolation(xStart,yStart,xEnd,yEnd,xValue):
    """ Takes in five arguments of type int.  xStart and yStart define corresponding posist\
        ions and values where the linear interpolation starts and xEnd and yEnd define correspondin\
        g posistions and values where the linear interpolation ends.  Returns the yValue that corre\
        sponds with the xValue parameter based on the defined line. """
    slope = (yEnd-yStart)/(xEnd-xStart)
    yValue = slope*(xValue-xStart)+yStart
    return yValue

def findChromNumb(file):
        """ Takes in a wig file and returns the number of chromosomes the file desc\
        ribes in integer form. """
        f = open(file,'r')
        text = f.read()
        chromNumb = text.count('span')
        return chromNumb
                    


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

def gapFunction2(location,dataDict,highestNumb,span):
    """Takes in a location (base#,chrom#) that isn't currently in the dataDict

    """
    chrom = location[1]
    #Finding lower bound
    lowerBound = location
    while True:
        if lowerBound in dataDict:
            break
        loc = lowerBound[0]
        loc+=-1
        lowerBound = (loc,lowerBound[1])
        if lowerBound[0] == 0:
            return dataDict

    #Finding upper bound
    upperBound = location
    while True:
        if upperBound in dataDict:
            break
        loc = upperBound[0]
        loc+=1
        upperBound = (loc,upperBound[1])
        if upperBound[0]>highestNumb:
            return dataDict

    #Calculating gap size
    gapSize = upperBound[0]-(lowerBound[0]+1)
    if gapSize<2*span:
        gap = 'small'
        #calculating value that should fill gap using linear Interpolation
        midx = (upperBound[0]+lowerBound[0]+1)/2.0
        lowerValue = float(dataDict[lowerBound])
        upperValue = float(dataDict[upperBound])
        i=0
        value = lowerValue
        while value == lowerValue:
            i+=1
            try: value = float(dataDict[(lowerBound[0]-i,chrom)])
            except: break
        xStart = lowerBound[0]-i+1
        xEnd = upperBound[0]
        midy = linInterpolation(xStart,lowerValue,xEnd,upperValue,midx)
        #add values to dataDict
        for x in range(gapSize):
            dataDict[(lowerBound[0]+1+x,chrom)]=midy
    return dataDict
            

def bedToWig(oldwig,tmpBed,newWig):

    bedFile = open(tmpBed,'r')
    bedLines = bedFile.readlines()
    span =findSpan(oldwig)
    
    #iterativing through old bed file and creating bedDict where {chrom#:(value,position)} and list of chromosomes.

    bedDict = {}
    chromList = []
    
    for l in bedLines:
        line = l.replace('\n','')
        dataArray = string.split(line, '\t')
        if int(dataArray[2])>int(dataArray[1]):
            start = int(dataArray[1])
            stop = int(dataArray[2])
        else:
            start = int(dataArray[2])
            stop = int(dataArray[1])
        chrom = dataArray[0]
        value = dataArray[4]
        
        if chrom in chromList:
            pass
        else:
#            print 'Obtaining data  for', chrom
            bedDict[chrom]=[]
            chromList.append(chrom)
        for location in range(start,stop):
            if location%span == 0:
                bedDict[chrom].append((value,location))

    #creating new file
    wigFile = open(newWig,'w')
    oldWigFile = open(oldWig, 'r')
    oldLines = oldWigFile.readlines()
    wigFile.write(oldLines[0])
    n_old = 0
    for chrom in chromList:
 #       print 'adding ', chrom, ' to new wig file'
        #writing the header line for each chromosome
        for n in range(n_old,len(oldLines)):
            if chrom in oldLines[n]:
                wigFile.write(oldLines[n])
                n_old = n
                break
        #writing the data for each chromosome
        sortedList = sorted(bedDict[chrom], key=operator.itemgetter(1))
        for n in range(len(sortedList)):
            wigFile.write(str(sortedList[n][1])+'\t'+sortedList[n][0]+'\n')
    return
     

if __name__ == '__main__':
    oldWig = sys.argv[1]
    tmpBed = sys.argv[2]
    newWig = sys.argv[3]
    bedToWig(oldWig,tmpBed,newWig)
