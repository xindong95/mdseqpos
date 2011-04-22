from Consts import NEG_INF
import array

def score2PeakInput(data):
    midpos_dict = {}
    score_dict = {}
    
    for seqId in data:
        midpos_dict[seqId] = array.array("I",[])
        m_append = midpos_dict[seqId].append
        score_dict[seqId] = array.array("d",[])
        s_append = score_dict[seqId].append
        for record in data[seqId]:
            (position, probeLength, ratios, numOfNearbyProbes, scoreValue, pValue) = record
            
            m_append(position+probeLength/2)
            if scoreValue:
                s_append(float(scoreValue))
            else:
                s_append(NEG_INF)
    
    return (midpos_dict,score_dict)

def makeNegPeakInput(peakInput, nullMean):
    score_dict = {}
    
    for seqId in peakInput[1]:
        score_dict[seqId] = array.array("d",[])
        
        pos_score_list = peakInput[1][seqId]
        append = score_dict[seqId].append
        for score in pos_score_list:
            if score == NEG_INF:
                negScore = NEG_INF
            else:
                negScore = 2 * nullMean - score
            
            append(negScore)
    
    return (peakInput[0],score_dict)

