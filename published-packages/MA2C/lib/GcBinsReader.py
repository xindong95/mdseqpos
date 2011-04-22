import csv
import math
from MA2C import CsvFilter
from MA2C import Consts

def _appendSignal(store, key, signal):
    
    if not store.has_key(key):
        store[key] = []
        
    store[key].append(math.log(signal, 2))
    

def _mergeBin(store, gcKeys, middle):
    # merge from lowest to middle
    if len(gcKeys) > 2:
        gcIndex = 0
        
        while True:
            gc = gcKeys[gcIndex]
            mergeToGc = gcKeys[gcIndex + 1]
            
            if mergeToGc > middle:
                break
            
            if len(store[gc]) < Consts.NORMALIZE_MINBIN:
                store[mergeToGc].extend(store[gc])
            
            gcIndex = gcIndex + 1
        
    # merge from highest to middle
    if len(gcKeys) > 2:
        gcIndex = len(gcKeys) - 1
        
        while True:
            gc = gcKeys[gcIndex]
            mergeToGc = gcKeys[gcIndex - 1]

            if mergeToGc < middle:
                break
            
            if len(store[gc]) < Consts.NORMALIZE_MINBIN:
                store[mergeToGc].extend(store[gc])
            
            gcIndex = gcIndex - 1
        

def read(tpmapFileName, rawFileName):
    ipGcBin = {}
    inputGcBin = {}
    
    rawReader = csv.DictReader(CsvFilter.Filter(rawFileName), delimiter = '\t')
    tpmapReader = csv.DictReader(CsvFilter.Filter(tpmapFileName), delimiter = '\t')
    
    for rawRow in rawReader:
        x = int(rawRow["X"])
        y = int(rawRow["Y"])
        ip = float(rawRow["IP"])
        input = float(rawRow["INPUT"])
        
        # match record in tpmap
        tpmapGc = None
        
        while True:
            tpmapRow = tpmapReader.next()
            
            tpmapX = int(tpmapRow["X"])
            tpmapY = int(tpmapRow["Y"])
            tpmapGc = int(tpmapRow["GC"])
            
            if (x, y) == (tpmapX, tpmapY):
                break
        
        _appendSignal(ipGcBin, tpmapGc, ip)
        _appendSignal(inputGcBin, tpmapGc, input)
        
    del tpmapReader
    del rawReader

    gcKeys = ipGcBin.keys()
    gcKeys.sort()
    
    lowest = gcKeys[0]
    highest = gcKeys[-1]
    middle = (lowest + highest) / 2.0
    
    _mergeBin(ipGcBin, gcKeys, middle)
    _mergeBin(inputGcBin, gcKeys, middle)
    
    return (ipGcBin, inputGcBin)