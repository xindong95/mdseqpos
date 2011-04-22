[data]
# Directories for bpmap and CEL files
# GenomeGrp: Genome Group in bpmap files. Human Chr21/22 and Encode arrays: (blank);
# Human Tiling 1.0 and 2.0 arrays: Hs;  Mouse arrays: Mm; Drosophila arrays: Dm
# bpmap and cel files should be in the same order and have the same keys (item before =)
# RepLib: full name of the the Repeat library file (optional)
# Group: separates all your cel files into Treatment (1) or Control (0) groups
#    (e.g. 111000 for 3ChIP and 3Input, 10 for 1ChIP 1Ipnut, etc)
# Paired: Specify 1 if the treatment and control experiments are strictly paired, otherwise leave it blank.
# If paired, the cel files should be well aligned, like: ChIP1 ChIP2 ChIP3 Ctrl1 Ctrl2 Ctrl3

BpmapFolder = /misc/iris/acct/weili/database/BPmap/Chr2122
CelFolder = /misc/iris/acct/weili/Tiling.Soft/Model.Tiling/lab/ER
GenomeGrp =
RepLib = /misc/iris/acct/weili/database/humanhg17_May2004/Annotation/Humanhg17Rep.lib
Group = 111000
Pair =

[bpmap]
1 = P1_CHIP_A.Anti-Sense.hs.NCBIv35.NR.bpmap
2 = P1_CHIP_B.Anti-Sense.hs.NCBIv35.NR.bpmap
3 = P1_CHIP_C.Anti-Sense.hs.NCBIv35.NR.bpmap

[cel]
1 = MCF_ER_A1.CEL  MCF_ER_A3.CEL  MCF_ER_A4.CEL MCF_INP_A1.CEL  MCF_INP_A3.CEL  MCF_INP_A4.CEL
2 = MCF_ER_B1.CEL  MCF_ER_B3.CEL  MCF_ER_B4.CEL MCF_INP_B1.CEL  MCF_INP_B3.CEL  MCF_INP_B4.CEL
3 = MCF_ER_C1.CEL  MCF_ER_C3.CEL  MCF_ER_C4.CEL MCF_INP_C1.CEL  MCF_INP_C3.CEL  MCF_INP_C4.CEL



[intensity analysis]
# BandWidth: 	the number of bases to extended from the position being analyzed.
# The result is that 2*Bandwidth probe positions are included in the signal and p-value analysis
# MaxGap: 	maximum gap between positive probes; All regions separated by < MaxGap will be mergered into one
# MinProbe: 	minimum number of probes for MAT score analysis

BandWidth = 	300
MaxGap = 	300
MinProbe  = 	10


[interval analysis]
# interval cutoff can be from either MAT score or pvalue or region False Discovery Rate (%) (FDR)
# if multiple cutoffs are specified, only one will be used as the real cutoff, others will be ignored.
# cutoff priority: FDR > Pvalue > Matscore
# default cutoff: Pvalue 1e-5
# Extend: extend the ChIP-enriched regions at both ends, default is 0.
# set it a negative number if you want to shrink the ChIP regions.

Matscore =
Pvalue = 1e-5
FDR =
Extend =

[output]
# the name of this TAG file will be used as the project name.
# Log: the screen output will be stored in a file (default False)

Log =








