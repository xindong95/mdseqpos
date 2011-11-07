#!/opt/bin/python
# Time-stamp: <2011-11-08 01:08:14 sunhf>
import sys
import os
import subprocess
from optparse import OptionParser
import os.path as Path
opt={}
run = subprocess.call
slim= lambda x: Path.splitext(Path.split(x)[1])[0]
tmp=lambda x:str(os.getpid())+"_."+slim(x)
'''
>>> tmp("/a/b/c.py)
'c'
'''
def check():
    with open(os.devnull,'w') as null:
        if os.waitpid(subprocess.Popen('bigWigAverageOverBed',stdout=null,stderr=null,shell=True).pid,0)[1] != 65280:
            print "Please add bigWigAverageOverBed to your Path"
            sys.exit(1)
        else:
            print "bugWigAverageOverBed checked successfully!"
def prepare_optparser():
    usage = "usage: %prog [options] <-r rfile> <-b bed file> <-w bigwig file>(>=2)"
    description = """Draw correlation plot for many bigwig files at regions by a bed file.

Method: It will calculate a value for each region defined in a bed
file based on each bigwig files. The method can be chosen from -m
option.
    """
    optparser = OptionParser(version="%prog 0.2",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-r","--rfile",dest="rfile",default="NA",
                         help="R output file. If not set, do not save R file.")
    optparser.add_option("-b","--bed",dest="bed",type="string",
                         help="the bed file you want to include in the calculation.")
    optparser.add_option("-w","--bw",dest="wig",type="string",action="append",
                         help="the bigwig file you want to include in the calculation. This option should be used for at least twice.",default=[])
    (options,args) = optparser.parse_args()
    global opt
    opt=options
    print opt
    if not options.wig or len(options.wig) < 2 or not options.rfile or not options.bed:
        optparser.print_help()
        sys.exit(1)
    # check the files
    if not os.path.isfile(options.bed):
        error("%s is not valid!" % options.bed)
        sys.exit(1)
    for f in options.wig:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)
    return options
def main():
    prepare_optparser()
    tmp_list,tmp_wig1_list,tmp_wig2_list=[],[],[]

    # convert the BED file into 4-lines style
    with open(opt.bed) as BED:
        for line in BED:
            tmp_list.append("\t".join(line.split("\t")[0:4])+"\n")

    # save the result into a temp bed file
    with open(tmp("bed"),'w') as TBED:
        TBED.writelines(tmp_list)

    # use bigWigAverageOverBed executive on bigwig files
    for a_wig in opt.wig:
        run(["bigWigAverageOverBed",a_wig,tmp("bed"),tmp(a_wig)])

    # extract the useful part from the output table
    cut_=lambda field:[line.split("\t")[field] for line in F]
    wig_average_profile_list=[]
    for a_wig in opt.wig:
        with open(tmp(a_wig)) as F:
            wig_average_profile_list.append(cut_(4))
            
    cj = lambda x:",".join(map(str,x))
    quoted_cj = lambda x:",".join(map(lambda x:"'"+str(x)+"'",x))    
    with open(tmp("r"),'w') as R:
        print >> R, "pdf('%s')"%tmp("pdf")
        for wig_num in range(len(wig_average_profile_list)):
            print >> R, "p%s <- c(%s)\n"%(wig_num,cj(wig_average_profile_list[wig_num]))

        print >> R, "c <- cbind(%s)"%(cj(["p"+str(i) for i in range(len(wig_average_profile_list))]))
        print >> R,'''
panel.plot <- function( x,y, ... ){
   par(new=TRUE)
   m <- cbind(x,y)
   plot(m,col=densCols(m),pch=20)
   lines(lowess(m[!is.na(m[,1])&!is.na(m[,2]),]),col="red")  
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(0, 1, 0, 1))
   r <- cor(x, y,use="complete.obs")
   txt <- format(round(r,2),width=5,nsmall=2)
   txt <- paste(prefix, txt, sep="")
   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
   text(0.5, 0.5, txt, cex = cex.cor)
}
pairs(c, lower.panel=panel.plot, upper.panel=panel.cor,labels=c(%s))
'''%(quoted_cj(map(slim,opt.wig)))
    print tmp("r")
    run(["Rscript",tmp("r")])
    print dir(opt)
    if opt.rfile!="NA":
        os.rename(tmp('r'),opt.rfile)
        os.rename(tmp('pdf'),opt.rfile+".pdf")
        print "\t= = The results are %s and %s = = "%(opt.rfile,opt.rfile+".pdf")
    else:
        print "\t= = The results are %s and %s = ="%(tmp('r'),tmp('pdf'))
        

def clean_():
    rm_ = lambda x:os.remove(x) if os.access(x,os.F_OK) else False
    rm_(tmp("bed"))
    print opt
    for a_wig in opt.wig:
        rm_(tmp(a_wig))

if __name__ == '__main__':
    try:
        check()
        main()
    except KeyboardInterrupt:
        clean_()
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
    except:

        raise
