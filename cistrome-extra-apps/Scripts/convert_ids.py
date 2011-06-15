#!/usr/bin/env python
# Time-stamp: <2011-06-15 16:48:46 Tao Liu>
"""
Convert gene ids through bioconductor.
"""
import subprocess
import os
import sys

rprog="""
sink(file=NULL,type='message')
suppressPackageStartupMessages(library(%s))

idfile  <- '%s'
outfile <- '%s'

genelist <- as.character(as.vector(read.table(idfile)[[1]]))

# conversion 
%s

converted <- converted[!is.na(converted)]

# write to output
write.table(converted,outfile,quote=F,col.names=F,row.names=F)
#

"""

def map_to_gene_universe ( dbkey ):
    """Mapping rule
    hg18 org.Hs.eg
    hg19 org.Hs.eg
    mm8	org.Mm.eg
    mm9	org.Mm.eg
    ce4	org.Ce.eg
    ce6	org.Ce.eg
    dm2	org.Dm.eg
    dm3	org.Dm.eg
    """
    if dbkey.startswith("hg"):
        # human
        return "org.Hs.eg"
    elif dbkey.startswith("mm"):
        # mouse
        return "org.Mm.eg"
    elif dbkey.startswith("ce"):
        # C. elegans
        return "org.Ce.eg"
    elif dbkey.startswith("dm"):
        # Fruitfly
        return "org.Dm.eg"
    else:
        raise Exception("Not a valid dbkey!")

def main():
    """
    convert_ids.py '$idfile.dbkey' '$idfile' '$conversion' '$outfile'
    """
    if len(sys.argv) < 5:
        sys.stderr.write("need 4 paras: %s <genome assembly> <id file> <conversion rule> <outfile>\n\n" % sys.argv[0])
        sys.stderr.write("It only supports mouse(mmN), human(hgN), C. elegans(ceN) or Drosophila Melanogaster(dmN)\n")
        sys.stderr.write("Rules should be either E2R, E2S, R2E, R2S, S2E, S2R where E stands for EntrezID, R stands for Refseq ID, and S for gene symbols.\n")
        sys.exit()

    gene_universe = map_to_gene_universe(sys.argv[1])
    input_file = sys.argv[2]
    output_file = sys.argv[4]
    rule = sys.argv[3]

    convert_script=""
    if rule == "E2R":
        convert_script += "x<-"+gene_universe+"REFSEQ\n"
        convert_script += "converted <- grep('NM',unlist(as.list(mget(genelist,x,ifnotfound=NA)),use.names=F),value=T)\n"
        convert_script += "converted <- converted[!is.na(converted)]\n"
    elif rule == "R2E":
        convert_script += "x<-"+gene_universe+"REFSEQ2EG\n"
        convert_script += "converted <- unlist(as.list(mget(genelist,x,ifnotfound=NA)),use.names=F)\n"
        convert_script += "converted <- converted[!is.na(converted)]\n"        
    elif rule == "E2S":
        convert_script += "x<-"+gene_universe+"SYMBOL\n"
        convert_script += "converted <- unlist(as.list(mget(genelist,x,ifnotfound=NA)),use.names=F)\n"
        convert_script += "converted <- converted[!is.na(converted)]\n"        
    elif rule == "S2E":
        convert_script += "x<-"+gene_universe+"SYMBOL2EG\n"
        convert_script += "converted <- unlist(as.list(mget(genelist,x,ifnotfound=NA)),use.names=F)\n"
        convert_script += "converted <- converted[!is.na(converted)]\n"        
    elif rule == "S2R":
        # first convert SYMBOL to ENTREZ
        convert_script += "x<-"+gene_universe+"SYMBOL2EG\n"
        convert_script += "tmpentrez <- unlist(as.list(mget(genelist,x,ifnotfound=NA)),use.names=F)\n"
        convert_script += "tmpentrez <- converted[!is.na(tmpentrez)]\n"        
        # then from ENTREZ to REFSEQ
        convert_script += "x<-"+gene_universe+"REFSEQ\n"
        convert_script += "converted <- grep('NM',unlist(as.list(mget(tmpentrez,x,ifnotfound=NA)),use.names=F),value=T)\n"
        convert_script += "converted <- converted[!is.na(converted)]\n"        
    elif rule == "R2S":
        # first convert REFSEQ to ENTREZ
        convert_script += "x<-"+gene_universe+"REFSEQ2EG\n"
        convert_script += "tmpentrez <- unlist(as.list(mget(genelist,x,ifnotfound=NA)),use.names=F)\n"
        convert_script += "tmpentrez <- converted[!is.na(tmpentrez)]\n"                
        # then from ENTREZ to SYMBOL
        convert_script += "x<-"+gene_universe+"SYMBOL\n"
        convert_script += "converted <- unlist(as.list(mget(tmpentrez,x,ifnotfound=NA)),use.names=F)\n"
        convert_script += "converted <- converted[!is.na(converted)]\n"                
    else:
        raise Exception("Unrecognized conversion %s" % rule)
    p = subprocess.Popen("R --vanilla", shell=True,executable="/bin/bash", stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    t = rprog % (gene_universe+".db",input_file,output_file,convert_script)
    (output,error) = p.communicate(t)
    if p.returncode != 0:
        sys.stderr.write("Error occurs during conversion! You may need to revise your parameters.\n")
        sys.stderr.write("Output: %s\n" % output)
        sys.stderr.write("Error: %s\n" % error)        

if __name__ == "__main__":
    main()

