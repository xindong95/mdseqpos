from MA2C.Functions import qnorm, pnorm

def calculate ( data , nullmean, nullstderr ):
    for seqId in data.keys():
        seqData = data[seqId]
        _dealSeqData(seqData,nullmean,nullstderr)
    return data


def _dealSeqData(seqData,nullmean, nullstderr):
    for i in xrange(len(seqData)):
        _dealOneRecord(seqData, i, nullmean, nullstderr)

def _dealOneRecord(seqData, i, nullmean, nullstderr):
    score = seqData[i][4]
    if score:
        r = pnorm((score - nullmean) / nullstderr, True);
        if (r == 0):
            seqData[i].append( 1e-300 )
        else:
            seqData[i].append( r )
    else:
        seqData[i].append( None)
