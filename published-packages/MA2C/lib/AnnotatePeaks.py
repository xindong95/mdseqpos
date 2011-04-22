from MA2C import Consts

def calculate(data, ppeaks, npeaks, FDR_table, bandwidth):
    
    ppeak_result = _annotate_peak_region(data, ppeaks, FDR_table, bandwidth)
    npeak_result = _annotate_peak_region(data, npeaks, FDR_table, bandwidth)
    
    return (ppeak_result, npeak_result)


def _annotate_peak_region(data, peaks, FDR_table, bandwidth):
    
    result = {}
    
    for seqId in peaks:
        result[seqId] = []
        append = result[seqId].append
        temp_peak = None        # [start_index, stop_index]
        for (start_index, stop_index, peak_score) in peaks[seqId]:
            
            if temp_peak:       # exists peak
                temp_peak_stop = data[seqId][temp_peak[1]][0]+data[seqId][temp_peak[1]][1]-1+bandwidth
                this_peak_start = data[seqId][start_index][0]-bandwidth
                if temp_peak_stop>this_peak_start: # merge!
                    temp_peak[1] = stop_index
                    continue
                else:
                    append( _close_annotate_peak (temp_peak[0],temp_peak[1],data[seqId], bandwidth, FDR_table) )
                    # close temp_peak
                    temp_peak = [start_index,stop_index]
            else:               # 
                temp_peak = [start_index,stop_index]
                continue
        if temp_peak:
            append( _close_annotate_peak (temp_peak[0],temp_peak[1],data[seqId],bandwidth, FDR_table) ) # for the last one
    return result

def _close_annotate_peak(start_index, stop_index,chrom_data, bandwidth, FDR_table):
    start_probe = chrom_data[start_index]
    stop_probe = chrom_data[stop_index]
            
    peak_start_position =  start_probe[0]-bandwidth
    peak_stop_position  =  stop_probe[0] + stop_probe[1] - 1+bandwidth
            
    (peak_summit_position, peak_pvalue, peak_score) = _find_summit_position(start_index, stop_index,chrom_data)
            
    peak_length = peak_stop_position - peak_start_position + 1
    peak_probe_num = stop_index - start_index + 1
    peak_summit_relative_pos = peak_summit_position-peak_start_position+1

    peak_FDR = _calc_FDR(peak_score, FDR_table)
            
    return (peak_start_position, peak_stop_position, 
            peak_summit_relative_pos, peak_length, 
            peak_probe_num, peak_score, peak_pvalue, peak_FDR)


def _find_summit_position(start_index, stop_index, chrom_data):
    
    summit_score = chrom_data[start_index][4]
    summit_index = start_index

    for i in xrange(start_index + 1, stop_index + 1):
        score = chrom_data[i][4]
        
        if score > summit_score:
            summit_index = i
            summit_score = score
    
    summit_probe = chrom_data[summit_index]
    
    summit_position = summit_probe[0] + summit_probe[1] / 2
    assert summit_score == summit_probe[4]

    summit_pvalue = summit_probe[5]
    
    return (summit_position, summit_pvalue, summit_score)


def _calc_FDR(score, FDR_table):
    result = Consts.UPPER_LIMIT_FDR
    
    for (fdr, score_value, pvalue, pos_len, neg_len) in FDR_table:
        if score_value <= score:
            if result > fdr:
                result = fdr

    return result
