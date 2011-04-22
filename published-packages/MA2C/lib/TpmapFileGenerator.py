from MA2C import PosFileReader
from MA2C import NdfFileReader
import os

def generate(tagFileInfo):
    posFileInfo = {}

    if not os.path.exists(tagFileInfo["output"]["dir"]):
        os.mkdir(tagFileInfo["output"]["dir"])
    fileName = tagFileInfo["output"]["tpmap"]

    # compare the modification time. If tpmap is newer than pos and ndf file, then skip building it.
    if os.path.exists(fileName)\
            and os.stat(tagFileInfo["sample"]["ndfFile"]).st_mtime < os.stat(fileName).st_mtime:
        if tagFileInfo["sample"].has_key("posFile"):
            if os.stat(tagFileInfo["sample"]["ndfFile"]).st_mtime < os.stat(fileName).st_mtime:
                return True
        else:
            return True

    if tagFileInfo["sample"].has_key("posFile"):
        posFileInfo = PosFileReader.read(tagFileInfo["sample"]["posFile"])
    
    ndfFileContent = NdfFileReader.read(tagFileInfo["sample"]["ndfFile"], posFileInfo)

    f = open(fileName, "w")
    
    try:
        f.write("#seq_group_name %s\n" % tagFileInfo["sample"]["designId"])
        f.write("#XYdimension\t%d\t%d\n" % (ndfFileContent["maxX"], ndfFileContent["maxY"]))
        # write the header line
        f.write("PROBE_SEQUENCE\tSTRAND\tSEQ_ID\tPOSITION\tX\tY\tGC\n")
        for row in ndfFileContent["tpmap"]:
            f.write("%s\t0\t%s\t%d\t%d\t%d" % row)
            # when writing tpmap, calculate GC counts for each probe 
            f.write("\t%d\n" % (row[0].count("C")+row[0].count("G")) )
                    
    finally:
        f.close()
    
    return True
