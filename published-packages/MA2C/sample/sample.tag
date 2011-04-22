[sample] 

# The DESIGN ID is unique per run. You can find this ID in your
# samplekey file from NimbleGen.

DESIGN_ID = 123

# NDF_FILE and POS_FILE: They are unique per run. Normally, you can
# find these files in DesignFiles directory. Valid path name either
# relative path or absolute path is accepted. You'd better keep these
# files consistent with DESIGN_ID (not strictly required).

NDF_FILE = 123.ndf
POS_FILE = 123.pos

# CHIP_ID can be found in your samplekey file from NimbleGen. A
# CHIP_ID means a single ChIP experiment with a pair of CYE5 and CYE3
# data. If multiple CHIP_IDs are assigned, MA2C will regard them as
# replicate experiments.

CHIP_ID  = 456 789

# IP_FILE and INPUT_FILE are the raw data for ChIP and input control
# data from NimbleGen, which can be found in your PairData
# directory. Sometimes, people applied dye swaps for their
# experiments, if so, cy5 will be the input and cy3 will be the IP in
# samplekey file from NimbleGen. Valid path name either relative path
# or absolute path is accepted. If multiple CHIP_IDs are assigned, the
# order in these two parameters MUST be consistent with the order for
# CHIP_ID parameter.

IP_FILE = 456_5_pair.txt 789_5_pair.txt
INPUT_FILE = 456_3_pair.txt 789_3_pair.txt

[normalization]
# Method: normalization method
# 		  Two options available: Robust, Simple

METHOD = Robust

# C     : a parameter to Robust method

C = 2

[peak detection]

# Method: criterion used for detecting ChIP-enriched regions. Can be
#         Pvalue or FDR or MA2C

METHOD = pvalue

# Threshold: cutoff score used in the above method.
#            For a 5% FDR, put 5; for 0.001 p-value, put 1e-3...

THRESHOLD = 1e-3

# Bandwidth (in basepairs): half-width of the sliding window

BANDWIDTH = 500

# MIN_PROBES: minimum number of probes required in the sliding window
#    centered at each probe; a probe having fewer probes than this
#    required number in its window will be ignored in the analysis.

MIN_PROBES = 5

# MAX_GAP: (in basepairs): maximum distance allowed for joining two
#    qualifying probes.

MAX_GAP = 250
