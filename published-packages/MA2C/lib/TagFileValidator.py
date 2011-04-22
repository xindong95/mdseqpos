import os.path

def _validateFileExists(filename):
    if not os.path.exists(filename):
        raise IOError("%s not found!" % filename)
    
def _validateFilesExisting(tagFileInfo):
    _validateFileExists(tagFileInfo["sample"]["ndfFile"])
    if tagFileInfo["sample"].has_key("posFile"):
        _validateFileExists(tagFileInfo["sample"]["posFile"])
        
    for chipId in tagFileInfo["sample"]["chips"]:
        _validateFileExists(tagFileInfo["sample"]["chips"][chipId]["ipFile"])
        _validateFileExists(tagFileInfo["sample"]["chips"][chipId]["inputFile"])

def validate(tagFileInfo):
    _validateFilesExisting(tagFileInfo)
