import unittest
from standardize_wig import *

class TestWigStandardization(unittest.TestCase):

    def testFindSpan(self):
        self.assertEquals(findSpan('/home/jacqui/data/truncatedWig1.wig'),50)

    def testFindChromNumb(self):
        self.assertEquals(findChromNumb('/home/jacqui/data/truncatedWig1.wig'),1)
    
    def testChromDict(self):
        datalines=['1\t.34\n','11\t.65\n','21\t.21\n','hello','31\t.91\n']
        span=2
        dataDict={1:'.34',2:'.34',11:'.65',12:'.65',21:'.21',22:'.21'}
        stopLine='hello'
        newDatalines=['31\t.91\n']
        lastNumb=21
        self.assertEquals(chromDict(datalines,span),[dataDict,stopLine,newDatalines,lastNumb])

    def testChromDict2(self):
        datalines=['1\t.34\n','11\t.65\n','21\t.21\n','hello','31\t.91\n']
        span=0
        dataDict={}
        stopLine='hello'
        newDatalines=['31\t.91\n']
        lastNumb=21
        self.assertEquals(chromDict(datalines,span),[dataDict,stopLine,newDatalines,lastNumb])

    def testChangeSpan(self):
        self.assertEquals(changeSpan('variableStep chrom=chrI span=30',13),'variableStep chrom=chrI span=13')

    def testStandardizeWig(self):
        pass
    
        
if __name__ == '__main__':
    unittest.main()
