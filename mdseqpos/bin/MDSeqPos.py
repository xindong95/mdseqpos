#! /usr/bin/env python

"""
SYNOPSIS: This tool searches for known and new (denovo) motifs that are in
regions specified by a bed file.
"""

import os
import sys
import optparse
import time
import random
import shutil
import json

from django.template import loader, Context
from django.core.management import setup_environ

from mdseqpos.motif import MotifList
from mdseqpos.chipregions import ChipRegions
import mdseqpos.settings as settings

setup_environ(settings)
#os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

USAGE = """USAGE: MDSeqPos.py BEDFILE GENOME
Arguments:
   BEDFILE - regions file
   GENOME  - assembly which the regions pertain to, e.g. 'hg18', 'mm9', etc.
             as defined in BUILD_DICT in lib/settings.py
   
Options:
   -d - flag to run denovo motif search (default: False)
   -m - comma separated list of known motifs dbs to use in the motif search,
        e.g. -m pbm.xml,transfac.xml
   -n - name of the output XML file which stores new motifs found during a
        denovo search, e.g. -n foo.xml (default: denovo.xml)
   -p - pvalue cutoff for the motif signficance, (default: 0.001)
   -s - name of species to filter the results with--if multiple species, comma
        separate them, e.g. hs,mm,dm
   -w - width of the region to be scanned for motifs; depends on resoution of
        assay, (default: 600 basepairs)
   -v - verbose: print out debug information (default: False)
   --maxmotif - the maximum number of motifs to report (defaut: 100)
"""
_DEBUG = False
_OUTPUT_DIR = "results"

#REMOVE _DEBUG as a parameter
def read_known_motifs(motif_dbs, _DEBUG):
    """Given a list of xml file names, this function tries to load the motifs
    in those databases
    """
    DATA_DIR = os.path.join(settings.DEPLOY_DIR, 'database')

    known_motifs = MotifList()
    for db in motif_dbs:
        if _DEBUG: print "loading (time): %s (%s)" % (db, time.time())
        tmp = MotifList()
        tmp.from_xml_file(os.path.join(DATA_DIR, db))
        known_motifs.extend(tmp)
        if _DEBUG: print "load Complete (time): %s (%s)" % (db, time.time())
    return known_motifs

def sample(p):
    """Given an array of probabilities, which sum to 1.0, randomly choose a 'bin',
    e.g. if p = [0.25, 0.75], sample returns 1 75% of the time, and 0 25%;
    NOTE: p in this program represents a row in the pssm, so its length is 4"""
    r = random.random()
    i = 0
    while r > p[i]:
        r -= p[i]
        i += 1
    return i

def pssm_to_fasta(pssm, fastafilename, n=1000):
    """Generate a fastafile that approximates the base frequences of the pssm"""
    base = ('A', 'C', 'G', 'T')
    fastafile = open(fastafilename, 'w')
    for seqnum in xrange(n):
        seqstr = "".join(base[sample(p)] for p in pssm)
        print >> fastafile, ">%d" % seqnum
        print >> fastafile, seqstr
    fastafile.close()

def make_logo(motif, width, height, img_dir):
    """Generate motif logo, place it in the img_dir and return its file name.
    """
    # create directory for holding motif logos if directory does not
    # already exist
    if not os.path.exists(img_dir):
        os.makedirs(img_dir)
    # generate temporary FASTA file which approximates the motif PSSM; input for seqlogo
    fasta_file_name = os.path.join(_OUTPUT_DIR, 'temp.fa')
    pssm_to_fasta(motif.pssm, fasta_file_name)
    logo_file_name = os.path.join(img_dir, motif.id)
    # run the seqlogo command
    command = "%s -f %s -o %s -w %d -h %d -F PNG -c -n -Y" % \
              (os.path.join(settings.DEPLOY_DIR, 'weblogo', 'seqlogo'),
               fasta_file_name, logo_file_name, width, height)
    os.system(command)
    # delete the temporary FASTA file
    os.remove(fasta_file_name)
                    
    # return the filename of the motif logo, strip out the _OUTPUT_DIR
    logo_file_name += '.png'
    return logo_file_name.replace(_OUTPUT_DIR+"/",'')

def save_to_html(motifList):
    """Saves the list of motifs as an html file named 'mdseqpos_out.html'.
    The motif logos associated with mdseqpos_out.html will be stored
    under _OUTPUT_DIR/img i.e. results/img.
    """
    #make results and results/img
    if not os.path.exists(_OUTPUT_DIR):
        os.makedirs(_OUTPUT_DIR)
    if not os.path.exists(os.path.join(_OUTPUT_DIR, 'img')):
        os.makedirs(os.path.join(_OUTPUT_DIR, 'img'))

    #make the motif logos, put them in results/img
    for (i, m) in enumerate(motifList):
        #NOTE: logoImg isn't part of the class definition, but we add it on to the instance
        #in hopes that the Motif.to_json will know how to serialize it.
        motifList[i].logoImg = make_logo(motifList[i], 192 / 20.0, 120 / 20.0,
                                         os.path.join(_OUTPUT_DIR, 'img'))

    #create the output file--
    mdseqpos_out_file = open(os.path.join(_OUTPUT_DIR,"mdseqpos_out.html"), 'w')
    motifList_json = motifList.to_json() 
    t = loader.get_template('mdseqpos_out.html')
    c = Context({'motifList': motifList_json}) 
    mdseqpos_out_file.write(t.render(c))

    #copy over supporting files--i.e. everything in the django/static dir
    for (path, dirs, files) in os.walk(os.path.join(settings.DEPLOY_DIR, 'django', 'static')):
        for f in files:
            shutil.copyfile(os.path.join(path, f), os.path.join(_OUTPUT_DIR, f))

    #END save_to_html
            
def main():
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-d', '--denovo', default=False, action="store_true",
                      help="flag to run denovo motif search (default: False)")
    parser.add_option('-m', '--known-motifs', default=None,
                      help="comma separated list of known motifs dbs to use \
                      in the motif search, e.g. -m pbm.xml,transfac.xml")
    parser.add_option('-n', '--new-motifs', default='denovo.xml',
                      help="name of the output XML file which stores new \
                      motifs found during adenovo search, e.g. -n foo.xml \
                      (default: denovo.xml)")
    parser.add_option('-p', '--pval', default=0.001,
                      help="pvalue cutoff for motif significance, \
                      (default: 0.001)")
    parser.add_option('-s', '--species-list', default=None, help="name of \
                      species to filter the results with--if multuple \
                      species, comma-separate them, e.g. hs,mm,dm")
    parser.add_option('-w', '--width', default=600,
                      help="width of the region to be scanned for motifs; \
                      depends on resoution of assay, (default: 600 basepairs)")
    parser.add_option('-v', '--verbose', default=False, action="store_true",
                      help="flag to print debugging info (default: False)")
    parser.add_option('--maxmotif', default=100,
                      help="maximum number of motifs to report,(default: 100)")
    
    #parse the command line options
    (opts, args) = parser.parse_args(sys.argv)
    if len(args) < 3: parser.error("Please specify a bedfile and/or a genome")
    bedfile_name = args[1]
    genome = args[2]
    _DEBUG = opts.verbose
    
    #READ in the regions that are specified in the BED file
    if _DEBUG: print "read regions start time: %s" % time.time()
    #HERE we should rely on a standard package to read in bed files; stub it
    chip_regions = ChipRegions(bedfile_name, genome)
    if _DEBUG: print "read regions end time: %s" % time.time()

    #LOAD the motifs (both known and denovo)
    known_motifs, new_motifs = None, None
    if opts.known_motifs:
        motif_dbs = [x.strip() for x in opts.known_motifs.split(',')]
        known_motifs = read_known_motifs(motif_dbs, _DEBUG)

    if opts.denovo:
        if _DEBUG: print "starting denovo search...(time: %s)" % time.time()
        new_motifs = chip_regions.mdmodule(width=int(opts.width))
        if _DEBUG: print "completed denovo search...(time: %s)" % time.time()
        new_motifs.save_to_xml(os.path.join(_OUTPUT_DIR, opts.new_motifs))
        
    #Combine both known and new motifs
    all_motifs = None
    if known_motifs and new_motifs:
        all_motifs = MotifList(known_motifs + new_motifs)
    elif known_motifs:
        all_motifs = known_motifs
    else:
        all_motifs = new_motifs

    #Run seqpos stats on all_motifs
    if _DEBUG: print "starting seqpos stats...(time: %s)" % time.time()
    for m in all_motifs: m.seqpos(chip_regions, width=int(opts.width))
    if _DEBUG: print "completed seqpos stats...(time: %s)" % time.time()

    #CULL the results to see only the relevant results, and output
    sig_motifs = all_motifs.cull(float(opts.pval), int(opts.maxmotif))
    
    #filter by species?
    if opts.species_list:
        species_list = opts.species_list.split(',')
        sig_motifs = sig_motifs.filterBySpecies(species_list)

    save_to_html(sig_motifs)
    
if __name__ == '__main__':
    main()
