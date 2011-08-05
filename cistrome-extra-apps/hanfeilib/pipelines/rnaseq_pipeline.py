#!/usr/bin/env python
import sys
import subprocess
_print=lambda x:sys.stdout.write("-"*10+'Input command:'+"-"*10+x+'\n')
_config=lambda x:lambda y:(subprocess.call(x+" "+y,shell=True),_print(x+" "+y))
_tools={"fastqc":_config("fastqc"),
        "tophat":_config("tophat"),
        "cufflinks":_config("cufflinks"),
        "cuffdiff":_config("cuffdiff")}
_p={"treat_fastq":[],
    "control_fastq":[],
    "tophat_out_treat":"tophat_out_treat",
    "tophat_out_control":"tophat_out_control",
    "refer_gtf":"",
    "refer_genome":"",
    "_index":{"hg19":"/mnt/Storage/data/Bowtie/hg19"},
    "_gtf":{"hg19":"/mnt/Storage/data/RefGene/hg19refGene.good.all.gtf"}
       }

def test():
    treat_add=_p["treat_fastq"].append
    control_add=_p["control_fastq"].append

    treat_add("liver10000.fastq")
    control_add("MAQC10000.fastq")
    _p["refer_genome"]=_p["_index"]["hg19"]
    _p["refer_gtf"]=_p["_gtf"]["hg19"]

    space_join=lambda L:' '.join(L)
    comma_join=lambda L:','.join(L)

    _tools["fastqc"](space_join(_p["treat_fastq"]))
    _tools["fastqc"](space_join(_p["control_fastq"]))

    map_back=lambda tophat_out,fastq:space_join(['-G',_p["refer_gtf"], '-o' ,_p[tophat_out],
                                 _p["refer_genome"],comma_join(_p[fastq])])

    _tools["tophat"](map_back("tophat_out_treat","treat_fastq"))
    _tools["tophat"](map_back("tophat_out_control","control_fastq"))

    bam_path=lambda tophat_out:_p[tophat_out]+"/accept_hits.bam"
    
    diff_gene=lambda treat_paths,control_paths:space_join([_p["refer_gtf"],
                                                           comma_join(treat_paths),
                                                           comma_join(control_paths)])
    _tools["cuffdiff"](diff_gene([bam_path("tophat_out_treat")],[bam_path("tophat_out_control")]))

test()

