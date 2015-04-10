# Time-stamp: <2011-07-22 13:45:10 Jian Ma>

''' script to read bigwig, only read meta information.'''

import struct
#uint2 = 'H'
#uint4 = 'I'
#uint8 = 'Q'

class BwIO:
    def __init__(self, filen = None):
        if filen:
            self.Read(filen)
    
    def PrintVars(self):
        arglist = [method for method in dir(self) if not callable(getattr(self, method))]
        for i in sorted(arglist):
            print "%s %s" %(i.ljust(30), getattr(self, i))
    
    def Read(self, filen):
        self.bwfh = open(filen, 'rb')
        bwfh = self.bwfh
        
        #bbiHeader, 64 bytes
        magic = bwfh.read(4)
        if magic == '&\xfc\x8f\x88':
            endianness = '<'
        elif magic == '\x88\x8f\xfc&':
            endianness = '>'
        else:
            raise IOError("The file is not in bigwig format")
        t = struct.unpack(endianness + 'HHQQQHHQQIQ', bwfh.read(60))
        k = {}
        k['version'], k['zoomLevels'], k['chromosomeTreeOffset'], k['fullDataOffset'], k['fullIndexOffset'], k['fieldCount'], k['definedFieldCount'], k['autoSqlOffset'], k['totalSummaryOffset'], k['uncompressBufSize'], k['reserved'] = t
        if k['version'] < 3:
            raise IOError("Bigwig files version <3 are not supported")
        self.bbiHeader = k
        
        #zoomHeaders, 24B, one for each zoomLevel
        self.zoomHeaders = []
        for i in range(self.bbiHeader['zoomLevels']):
            t = struct.unpack("IIQQ", bwfh.read(24))
            self.zoomHeaders.append({
                'reductionLevel': t[0],
                'reserved': t[1],
                'dataOffset': t[2],
                'indexOffset': t[3],
                })
        
        #totalSummary, 40B
        t = struct.unpack("Qdddd", bwfh.read(40))
        k = {}
        k['basesCovered'], k['minVal'], k['maxVal'], k['sumData'], k['sumSquares'] = t
        self.totalSummary = k

        #chromosomeTree, 
        bwfh.seek(self.bbiHeader['chromosomeTreeOffset'])
        magic = bwfh.read(4)
        if magic == '\x91\x8c\xcax':
            endianness = '<'
        elif magic == 'x\xca\x8c\x91':
            endianness = '>'
        else:
            raise ValueError("Wrong magic for this bigwig data file")
        t = struct.unpack(endianness + "IIIQQ", bwfh.read(28))
        k={}
        k['blockSize'], k['keySize'], k['valSize'], k['itemCount'], k['reserved'] = t
        
        t = struct.unpack(endianness + "bbH", bwfh.read(4))
        k['isLeaf'], k['reserved'], k['count'] = t
        self.chromosomeTree = k

        self.chromosomeTree['nodes'] = []
        for i in range(self.chromosomeTree['itemCount']):
            if self.chromosomeTree['isLeaf'] == 1:
                t = struct.unpack(endianness + "%dsII" %self.chromosomeTree['keySize'], bwfh.read(self.chromosomeTree['keySize'] + 8))
                self.chromosomeTree['nodes'].append({
                    'key': t[0].strip('\x00'),
                    'chromId': t[1],
                    'chromSize': t[2],
                    })
            elif self.chromosomeTree['isLeaf'] == 0:
                t = struct.unpack(endianness + "%dsQ" %self.chromosomeTree['keySize'], bwfh.read(self.chromosomeTree['keySize'] + 8))
                self.chromosomeTree['nodes'].append({
                    'key': t[0].strip('\x00'),
                    'childOffset': t[1],
                    })
        
        #dataCount
        #self.dataCount = struct.unpack("I", bwfh.read(4))[0]
        
'''  
        #data
        self.data = {}
        t = struct.unpack("IIIIIBBH", bwfh.read(24))
        self.data['chromId'], self.data['chromStart'], self.data['chromEnd'], self.data['itemStep'], self.data['itemSpan'], self.data['type'], self.data['reserved'], self.data['itemCount'] = t

        #index
        self.index = {}
        t = struct.unpack("IIQIIIIQII", bwfh.read(48))
        self.index['magic'], self.index['blockSize'], self.index['itemCount'], self.index['startChromIx'], self.index['startBase'], self.index['endChromIx'], self.index['endBase'], self.index['endFileOffset'], self.index['itemsPerSlot'], self.index['Reserved'] = t
        self.index['node'] = []
        print self.index['magic']
        print self.index['itemCount']
        for i in range(self.index['itemCount']):
            t = struct.unpack("bbH", bwfh.read(4))
            self.index['node'].append({
            'isLeaf': t[0],
            'reserved': t[1],
            'count': t[2],
            })
            if t[0] == 1:
                t2 = struct.unpack("IIIIQQ", bwfh.read(32))
                self.index['node'][-1]['startChromIx'], self.index['node'][-1]['startBase'], 
                self.index['node'][-1]['endChromIx'], self.index['node'][-1]['endBase'], 
                self.index['node'][-1]['dataOffset'], self.index['node'][-1]['dataSize'] = t2
            elif t[0] == 0:
                t2 = struct.unpack("IIIIQ", bwfh.read(24))
                self.index['node'][-1]['startChromIx'], self.index['node'][-1]['startBase'], 
                self.index['node'][-1]['endChromIx'], self.index['node'][-1]['endBase'], 
                self.index['node'][-1]['dataOffset'] = t2
            '''
            



        
'''
reload br
s=br.BwIO()


'''

