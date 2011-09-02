#!/usr/bin/python
import os
import sys
import re
import logging
import subprocess
import warnings
from optparse import OptionParser
import exonCEAS.inout as inout
import exonCEAS.exonplot as exonplot

logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

error   = logging.critical		# function alias
warn    = logging.warning
info    = logging.info

_print=lambda x:sys.stdout.write(x+'\n')
_p_stdout=lambda command:subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
run_cmd=lambda command:(subprocess.call(command,shell=True),_print(command))


    

def main():


    # read the options and validate them
    options=opt_validate(prepare_optparser())
    # read the gene annotation table
    info("# read the gene table...")
    GeneT = inout.GeneTable()
    GeneT.read(Host = None, User= None, Db=options.gdb, annotation='GeneTable', columns=('name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','exonStarts', 'exonEnds', 'name2'))
    GeneT.sort()

    chroms_GeneT=[i for i in GeneT.get_chroms() if not re.search('_[A-Za-z0-9]*',i)]
    info("# group gene and calculate profile")
    group_genes_and_plot_exons(options,GeneT,chroms_GeneT)

    try:
        p = subprocess.Popen("R" + " --vanilla < %s"  %(options.name+'_CI.R'), shell=True)
        p.wait()
        info ('#... Great! See %s for the graphical results of CEAS!' %(options.name+'_CI.pdf'))
    except:
        info ('#... Oops! Run %s using R for the graphical results of CEAS! CEAS could not run R directly.' %(options.name+'.R'))
    	

def read_wiggle(wiggle_path):

    wiggle_file=open(wiggle_path,"r")
    wiggle_content=wiggle_file.read()
    wiggle_chroms_re=re.compile(r'chrom=(\S+)\s')

    wiggle_chroms=wiggle_chroms_re.findall(wiggle_content)

    wiggle_splited=wiggle_content.split("chrom=")[1:]
    wiggle_profile_re=re.compile(r'(\d+)\t([\-]*\d+[\.]*\d*)')

    wiggle_dictionary={}

    for (one_section,one_chrom) in zip(wiggle_splited,wiggle_chroms):
        info( "reading wiggle file (%s)'s %s"%(wiggle_path,one_chrom))
        profile_list=wiggle_profile_re.findall(one_section)
        transposed_section=zip(*((int(i),float(j)) for (i,j) in profile_list))
        wiggle_dictionary[one_chrom]=transposed_section
    wiggle_file.close()
    return wiggle_dictionary


def group_genes_and_plot_exons(options,GeneT,chroms_GeneT):
                             
    span=max(options.ex_ispan,options.ex_espan)

    label=lambda x,y:options.name+{"peak":"_withpeak","flat":"_withoutpeak"}[x]+{"start":"_es","end":"_ed"}[y]
    path_bed=lambda x,y:label(x,y)+".bed"
    path_dump=lambda x,y:label(x,y)+"_dump.txt"
    

    write=lambda x:open(x,'w')
    
    if options.bed:
        info("Bed file(result of peak calling) is input.Auto group by bed peaks")
        _bed_read=lambda x:[i.strip().split() for i in open(x,"rU") if not i.startswith('#')]
        bed_list=_bed_read(options.bed)

        def read_p(command):
            p=_p_stdout(command)
            p.wait()
            stdout=p.stdout.read()
        
        info_BW=read_p("bigWigInfo "+options.wig+" -chroms")
        chroms_BW=re.compile("\r\n\t(chr\w+)").findall(info_BW)
        
        intersect=lambda x,y:set(x)&set(y)
        for chr in intersect(chroms_GeneT,chroms_BW):
            info("#getting genes on %s"%chr)
            ixs1,ixs2=exonplot.get_gene_indicies_by_bedpeaks(GeneT,bed_list,chr)

            list1=exonplot.get_exons_byindex(GeneT, ixs1, chr, options)
            list2=exonplot.get_exons_byindex(GeneT, ixs2, chr, options)
            if (list1 and list2) !=False:
                exonplot.paired_bed_make(list1, write(path('peak','start')),write(path('peak','end')) , chr)
                exonplot.paired_bed_make(list2, write(path('flat','start')),write(path('flat','end')) , chr)

        _prefix="siteproBW --dir --dump --span="+span+" --pf-res="+options.pf_res+" -w "+options.wig
        _suffix=lambda x,y:" -l "+label(x,y)+" -b "+path(x,y)

        run_cmd(_prefix+_suffix('peak','start')+_suffix('peak','end'))
        run_cmd(_prefix+_suffix('flat','start')+_suffix('flat','end'))
        
    bin=lambda span:span / options.pf_res
    Rscript_text=exonplot.make_Rscript_with_CI(options.gn_names,options.name,bin(options.ex_ispan), bin(options.ex_espan),options.pf_res,span)
    #use 'exonplot' to draw all groups in one page and calculate confidence interval        

    with open(options.name+'_CI.R','w') as r:
        r.write(Rscript_text)
        r.write('dev.off()')

    if not options.dump:
        
        rm=lambda:x:map(os.remove,x)
        rm([a(b,c) for a in (path_dump,path_bed) for b in ("peak","flat") for c in ("start","end")])

def opt_validate (optparser):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    # if gdb not given, print help, either BED or WIG must be given 
    if not options.gdb and not options.bed and not options.wig:
        optparser.print_help()
        sys.exit(1)
    elif not options.gdb:
        error('A gene table file must be given through -g (--gt).')
        sys.exit(1)
    elif options.gdb and not options.bed and not options.wig:
        error('Either a BED file or a WIG file must be given.')
        sys.exit(1)
   
    ##
    # check what inputs are given and determine which modules will operate
    ##
    
    #
    # check gene annotation table database
    # 
    HAVELOCALGDB = os.path.isfile(options.gdb)
    if not HAVELOCALGDB:
        error("No such gene table file as '%s'" %options.gdb)
        sys.exit(1)
    else:
        options.gdbtype = 'localdb'
        options.Host = None
        options.User = None
    
    #
    #check the ChIP bed file
    #
    if options.bed:
        HAVEBED = os.path.isfile(options.bed)
        if not HAVEBED:
            error("Check -b (--bed). No such bed file as '%s'" %options.bed)
            sys.exit(1)
        if os.path.getsize(options.bed) > 5000000:
            warnings.warn("ChIP bed file size may be too large to run CEAS with. Make sure it is a 'peak' file!")
            #error("ChIP bed file size is too big to handle! The file size is limmited to 5M bytes.")
            #sys.exit(1)
    else: HAVEBED = False
    
    #
    # check the wig file
    # 
    if options.wig:
        HAVEWIG=os.path.isfile(options.wig)
        if not HAVEWIG:
            error("Check -w (--wig). No such wig file as '%s'" %options.wig)
            sys.exit(1)
    else: HAVEWIG=False
        
    # get the experiment name
    #
    # if options.name is not given, BED and WIG file names will be used in order
    if not options.name:
        if HAVEBED:
            options.name=os.path.split(options.bed)[-1].rsplit('.bed',2)[0]
        elif HAVEWIG:
            options.name=os.path.split(options.wig)[-1].rsplit('.wig',2)[0]
    
    # profiling resolution
    return options
# ------------------------------------
# functions
# ------------------------------------
  
def prepare_optparser ():
    """Prepare optparser object. New options will be added in this
    function first.
    
    """
    
    usage = "usage: %prog < input files > [options]"
    description = "Exon CEAS (Cis-regulatory Element Annotation System for Exon)"
    
    optparser = OptionParser(version="%prog -- 0.9.9.8 beta (package version 1.0.2)",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-b","--bed",dest="bed",type="string",
                         help="BED file of ChIP regions.")
    optparser.add_option("-w","--wig",dest="wig",type="string",
                         help="WIG file for either wig profiling or genome background annotation. WARNING: --bg flag must be set for genome background re-annotation.")
    optparser.add_option("-g","--gt",dest="gdb",type="string",
                         help="Gene annotation table (eg, a refGene table in sqlite3 db format provided through the CEAS web, http://liulab.dfci.harvard.edu/CEAS/download.html).")
    optparser.add_option("--name",dest="name",\
                         help="Experiment name. This will be used to name the output files. If an experiment name is not given, the stem of the input BED file name will be used instead (eg, if 'peaks.bed', 'peaks' will be used as a name.)")
    optparser.add_option("--pf-res", dest="pf_res", type="int",\
                          help="Wig profiling resolution, DEFAULT: 50bp. WARNING: Value smaller than the wig interval (resolution) may cause aliasing error.", default=50) 
    optparser.add_option("--utr", action="store_true", dest="exon_utr",\
                         help="Whether to select the first and the last exon(next to 5' utr or 3' utr).If set,these exons will be included.",default=False)
    optparser.add_option("--espan", dest="ex_espan", type="int",\
                         help="Span from exon boundaries to exon region,DEFAULT=500bp", default=300)
    optparser.add_option("--ispan", dest="ex_ispan", type="int",\
                         help="Span from exon boundaries to ,DEFAULT=500bp", default=300)

    optparser.add_option("--dump", action="store_true", dest="dump",\
                     help="Whether to save the raw profiles of near exon boundary. The file names have a suffix of XXX, and YYY after the name.",default=False)
    return optparser
                    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        warn("User interrupts me! ;-) See you!")
        sys.exit(0)
