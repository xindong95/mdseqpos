# Time-stamp: <2008-06-03 17:03:50 Tao Liu>

"""Module Description: Read PairData and tpmap to write raw data.

"""

# ------------------------------------
# python modules
# ------------------------------------
import csv
import os
from MA2C import CsvFilter
from MA2C import ChipDataReader

# ------------------------------------
# Misc functions
# ------------------------------------


def _readTpmapFile(tpmapFileName):
    result = []
    
    tpmapReader = csv.DictReader(CsvFilter.Filter(tpmapFileName), delimiter = '\t')
    
    for row in tpmapReader:
        x = int(row["X"])
        y = int(row["Y"])
        
        result.append((x, y))
    
    del tpmapReader
    
    return result


def generate(tagFileInfo):

    result = []
    
    tpmapCache = None
    
    chips = tagFileInfo["sample"]["chips"]
    tpmapFileName = tagFileInfo["output"]["tpmap"]

    if not os.path.exists(tagFileInfo["output"]["dir"]):
        os.mkdir(tagFileInfo["output"]["dir"])

    for chipId in chips:
        chipRawFileName = tagFileInfo["output"]["chips"][chipId]["raw"]

        # compare modification time, if chipRawFile is newer than tpmap, ip/input file, then skip this chipId
        if os.path.exists(chipRawFileName)\
                and os.stat(tpmapFileName).st_mtime < os.stat(chipRawFileName).st_mtime\
                and os.stat(chips[chipId]["ipFile"]).st_mtime < os.stat(chipRawFileName).st_mtime\
                and os.stat(chips[chipId]["inputFile"]).st_mtime < os.stat(chipRawFileName).st_mtime:
            continue
        
        if tpmapCache == None:
            tpmapCache = _readTpmapFile(tpmapFileName)
            
        chipData = ChipDataReader.read(chips[chipId])
        
        f = open(chipRawFileName, "w")
        try:
            # header lines
            f.write("#DESIGN_ID\t%s\n" % tagFileInfo["sample"]["designId"])
            f.write("X\tY\tIP\tINPUT\n")
            
            for (x, y) in tpmapCache:
                key = (x, y)
                
                # must tolerate some missing data cases. just skip the
                # probe if it is missing in either ip or input.
                if not chipData["ip"].has_key(key) or not chipData["input"].has_key(key):
                    # give a warning and continue
                    continue

                ip = chipData["ip"][key]
                input = chipData["input"][key]
                f.write("%d\t%d\t%s\t%s\t\n" % (x, y, ip, input))
            
        finally:
            f.close()

    return True
