import csv
from MA2C import CsvFilter

def read(fileName):
    reader = csv.DictReader(CsvFilter.Filter(fileName), delimiter = '\t')
    
    result = {}
    
    for row in reader:
        seqId = row.setdefault("SEQ_ID", None)
        if seqId == None:
            raise IOError("SEQ_ID column does not exist in %s!" % fileName)
        
        chromosome = row.setdefault("CHROMOSOME", None)
        if chromosome == None:
            raise IOError("CHROMOSOME column does not exist in %s!" % fileName)
        
        result[seqId] = chromosome
    
    return result
        
