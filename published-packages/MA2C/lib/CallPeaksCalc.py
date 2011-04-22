# Time-stamp: <2008-05-22 17:26:10 Tao Liu>

"""Module Description: General Module to call peaks for microarry
data.

"""

# ------------------------------------
# python modules
# ------------------------------------
#from time import time

# ------------------------------------
# Misc functions
# ------------------------------------

def calculate(peak_input, threshold, max_gap):
    """General peak calling function.

    Arg: 

    1. peak_input is a tuple as (middle_point_dict,score_dict).
       
       middle_point_dict is a dict with seq_id as keys and lists of
       middle points for probes ( sorted by coordinates) as values.

       score_dict is a dict with seq_id as keys and lists of scores
       for probes ( sorted by coordinates) as values.

    2. threshold is a minimum score for peak.

    3. max_gap is the maximum gap to merge nearby probe to current
    peak.

    Return:

    peak_index_dict is a dictionary with seq_id as keys and lists of
    tuple (start_index,stop_index,MA2Cscore) as value.

    """
    result = {}
    #ta = time()
    midpos_dict = peak_input[0]
    score_dict = peak_input[1]
    
    for seq_id in midpos_dict:
        result[seq_id] = []
        append = result[seq_id].append
        chrom_midpos = midpos_dict[seq_id]
        chrom_score  = score_dict[seq_id]
        
        i = -1
        tpeak = None            # temporarily peak information holder
        for score_value in chrom_score:
            i = i+1

            if score_value <= threshold :
                if tpeak:
                    append(tpeak)
                    tpeak = None
            else:
                #check if tpeak exists
                if tpeak: # old peak
                    if chrom_midpos[i]-chrom_midpos[i-1] <= max_gap:
                        #extend
                        tpeak[1] = i
                        if score_value>tpeak[2]:
                            tpeak[2] = score_value
                    else:
                        append(tpeak)
                        tpeak = None
                else:        # new peak
                    tpeak=[i,i,score_value]
        if tpeak:               # last one
            append(tpeak)

    #print time()-ta
    return result

