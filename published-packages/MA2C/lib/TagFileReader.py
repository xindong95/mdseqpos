# Time-stamp: <2009-06-19 16:51:20 Tao Liu>
"""Module Description: Parse the tag file, store parameters, input and
output filenames

@status:  experimental
@version: $Revision$
"""


import ConfigParser
import string
from os.path import join as pjoin

class ConfigMissingException(Exception):
    pass

class IpNumberException(Exception):
    pass

class InputNumberException(Exception):
    pass

class WrongConfigValueException(Exception):
    pass

def read(tagFile):
    config = {}
    
    cp = ConfigParser.ConfigParser()

    raw_tag_name  = tagFile.rsplit('.',1)[0]
    
    if len(cp.read(tagFile)) == 0:
        raise IOError("%s not found!" % tagFile)
    
    for sec in cp.sections():
        secName = string.lower(sec)
        for opt in cp.options(sec):
            optName = string.lower(opt)
            config[secName + "." + optName] = string.strip(cp.get(sec, opt))
    
    result = {}
    
    result["tagFile"] = tagFile
    result["rawName"] = raw_tag_name
    result["sample"] = {}       # for input
    result["output"] = {"dir":"MA2C_Output"} # for output; the
                                             # directory name for
                                             # outpu

    # sample design id
    sampleDesignId = config.setdefault("sample.design_id", None)
    if sampleDesignId == None:
        raise ConfigMissingException("DESIGN_ID in SAMPLE section is missing!")
    
    result["sample"]["designId"] = sampleDesignId
    result["output"]["tpmap"] = pjoin("MA2C_Output","MA2C_%s.tpmap" % sampleDesignId) # tpmap filename

    # sample ndf file
    sampleNdfFile = config.setdefault("sample.ndf_file", None)
    if sampleNdfFile == None:
        raise ConfigMissingException("NDF_FILE in SAMPLE section is missing!")
    
    result["sample"]["ndfFile"] = sampleNdfFile
    
    # sample pos file
    samplePosFile = config.setdefault("sample.pos_file", None)
    if samplePosFile != None:
        result["sample"]["posFile"] = samplePosFile
    
    # sample chips, ipFiles & inputFiles  
    chipIds = config.setdefault("sample.chip_id", None)
    if chipIds == None:
        raise ConfigMissingException("CHIP_ID in SAMPLE section is missing!")
    
    chipIdList = chipIds.split()

    result["output"]["FDRtable"] = pjoin("MA2C_Output",raw_tag_name+"_FDRtable.txt")
    result["output"]["PeakXLS"] = pjoin("MA2C_Output",raw_tag_name+"_peaks.xls")
    result["output"]["NegPeakXLS"] = pjoin("MA2C_Output",raw_tag_name+"_negative_peaks.xls")
    result["output"]["PeakBED"] = pjoin("MA2C_Output",raw_tag_name+"_peaks.bed")
    result["output"]["wiggle"] = pjoin("MA2C_Output",raw_tag_name+"_MA2Cscore.wig.gz")

    
    ipFiles = config.setdefault("sample.ip_file", None)
    if ipFiles == None:
        raise ConfigMissingException("IP_FILE in SAMPLE section is missing!")
    
    ipFileList = ipFiles.split()
    
    inputFiles = config.setdefault("sample.input_file", None)
    if inputFiles == None:
        raise ConfigMissingException("INPUT_FILE in SAMPLE section is missing!")
    
    inputFileList = inputFiles.split()
    
    if (len(chipIdList) != len(ipFileList)):
        raise IpNumberException("The number of IP_FILE in SAMPLE section does not match the number of CHIP_ID!")
    
    if (len(chipIdList) != len(inputFileList)):
        raise InputNumberException("The number of INPUT_FILE in SAMPLE section does not match the number of CHIP_ID!")
    
    # put chips in the result
    result["sample"]["chips"] = {}
    result["output"]["chips"] = {}
    
    for i in xrange(len(chipIdList)):
        chipId = chipIdList[i]
        result["sample"]["chips"][chipId] = {}
        result["output"]["chips"][chipId] = {}
        result["sample"]["chips"][chipId]["ipFile"] = ipFileList[i]
        result["sample"]["chips"][chipId]["inputFile"] = inputFileList[i]
        result["output"]["chips"][chipId]["raw"] = pjoin("MA2C_Output","MA2C_%s_raw.txt" % chipId)
        result["output"]["chips"][chipId]["normalized"] = pjoin("MA2C_Output","MA2C_%s_normalized.txt" % chipId)
        
    result["normalization"] = {}
    
    method = config.setdefault("normalization.method", None)
    if method == None:
        raise ConfigMissingException("METHOD in NORMALIZATION section is missing!")
    method = string.lower(method)
    
    if method != "simple" and method != "robust":
        raise WrongConfigValueException("METHOD in NORMALIZATION section must be Simple or Robust!")
    
    result["normalization"]["method"] = method
    
    if (method == "robust"):
        c = config.setdefault("normalization.c", None)
        if c == None:
            raise ConfigMissingException("C in NORMALIZATION section is missing!")
        result["normalization"]["c"] = int(c)
    
    # peak detection
    result["peakDetection"] = {}
    
    method = config.setdefault("peak detection.method", None)
    if method == None:
        raise ConfigMissingException("METHOD in PEAK DETECTION section is missing!")
    
    result["peakDetection"]["method"] = string.lower(method)
    
    threshold = config.setdefault("peak detection.threshold", None)
    if threshold == None:
        raise ConfigMissingException("THRESHOLD in PEAK DETECTION section is missing!")
    
    result["peakDetection"]["threshold"] = float(threshold)
    
    bandwidth = config.setdefault("peak detection.bandwidth", None)
    if bandwidth == None:
        raise ConfigMissingException("BANDWIDTH in PEAK DETECTION section is missing!")
    
    result["peakDetection"]["bandwidth"] = int(bandwidth)
    
    minProbes = config.setdefault("peak detection.min_probes", None)
    if minProbes == None:
        raise ConfigMissingException("MIN_PROBES in PEAK DETECTION section is missing!")
    
    result["peakDetection"]["minProbes"] = int(minProbes)
    
    maxGap = config.setdefault("peak detection.max_gap", None)
    if maxGap == None:
        raise ConfigMissingException("MAX_GAP in PEAK DETECTION section is missing!")
    
    result["peakDetection"]["maxGap"] = int(maxGap)
    
    return result
