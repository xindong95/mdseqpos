#! /usr/bin/env python

"""A command-line script for running MDmodule and SeqPos."""

import os
import sys
import shutil
import time
import mdseqpos

from optparse import OptionParser
from mdseqpos.chipregions import ChipRegions
from mdseqpos.motif import MotifList, MotifTree

import mdseqpos.settings as settings

USAGE = """Usage: MDSeqPos.py [options] BEDFILE GENOME

Arguments:
  BEDFILE              name of an existing BED file containing ChIP regions,
                       to be scanned for motifs
  GENOME               abbreviation identifying the genome build that is
                       used in the BED file (must be 'hg18', 'mm3', or ... )"""

#DEPLOY_DIR = mdseqpos.__path__[0]
DATA_DIR = os.path.join(settings.DEPLOY_DIR, 'database')
STATIC_IMG_DIR = os.path.join(settings.DEPLOY_DIR, 'img')

#KNOWN_MOTIFS_FILE = 'known.xml'

RESULTS_DIR = 'results/'
NEW_MOTIFS_FILE = 'denovo.xml'
OUTPUT_FILE = 'mdseqpos_out.html'

if __name__ == '__main__':
    # parse command line arguments
    # 
    # Required Arguments:
    # 1. chip_regions_file_name -- BED file of ChIP regions
    # 2. genome -- genome identifier (e.g. 'hg18')
    # 
    # Optional Arguments:
    # 1. known_motifs_file_name -- input XML file where known motifs are stored
    # 2. new_motifs_file_name -- output XML file to store new motifs
    # 3. html_file_name -- name of output HTML file
    # 
    parser = OptionParser(usage=USAGE)
    parser.add_option('-m', '--known-motifs', default=None, help="name of input XML file containing known motifs, to be scanned against ChIP regions--if multuple xml files, comma-separate them, e.g. pbm.xml,transfac.xml")
    parser.add_option('-M', '--known-motifs-directory', default=DATA_DIR, help="directory of input XML file containing known motifs")
    parser.add_option('-n', '--new-motifs', default=NEW_MOTIFS_FILE, help="name of output XML file for storing new motifs discovered by de novo scan of ChIP regions")
    parser.add_option('-N', '--new-motifs-directory', default=RESULTS_DIR, help="directory of output XML file for storing new motifs")
    parser.add_option('-o', '--output', default=OUTPUT_FILE, help="name of output HTML file for display of MDSeqPos results")
    parser.add_option('-O', '--output-directory', default=RESULTS_DIR, help="directory of output HTML file")
    parser.add_option('-I', '--img-output-directory', default=os.path.join(RESULTS_DIR, "img"), help="directory of output motif logo img file")
    parser.add_option('-d', '--denovo', action="store_true", help="Flag to run denovo motif search (default: False)")
    parser.add_option('-p', '--pval', default=0.001, help="p-value cutoff for motif significance, 0.001 default") 
    parser.add_option('-w', '--width', default=600, help="width of region to be scanned for motif depends on resolution of assay, 600bp default") 
    parser.add_option('--maxmotif', default=0, help="maximum number of motifs to report. Default: 0, means output all motifs fit pval cut off.")
    parser.add_option('-s', '--species-list', default=None, help="name of species to filter the results with--if multuple species, comma-separate them, e.g. hs,mm,dm")

    parser.set_defaults(denovo=False)
    (opts, args) = parser.parse_args(sys.argv)
    if len(args) < 3:
        parser.error("must provide 2 required arguments specifying BEDFILE and GENOME")
    # required arguments
    chip_regions_file_name = args[1]
    genome = args[2]

    width = int(opts.width)
    pvalcutoff = float(opts.pval)
    maxmotif = int(opts.maxmotif)

    print 'start time', time.ctime()
    # retrieve ChIP regions from BED file
    chip_regions = ChipRegions(chip_regions_file_name, genome)
    print 'bed read time', time.ctime()

    known_motifs = None
    # handle/load known motifs
    if opts.known_motifs is not None:
        # known motifs directory and file name
        #if opts.known_motifs_directory[-1] != '/':
        #    opts.known_motifs_directory += '/'
        if not os.path.exists(opts.known_motifs_directory):
            os.makedirs(opts.known_motifs_directory)
        #Motif dbs are comma separated
        known_motifs_file_names = [os.path.join(opts.known_motifs_directory,
                                                name)
                                   for name in opts.known_motifs.split(",")]
        
        # retrieve known motifs from known motifs XML file
        known_motifs = MotifList()
        for motif_db in known_motifs_file_names:
            print "LOADING: %s" % motif_db
            tmp = MotifList()
            tmp.from_xml_file(motif_db)
            known_motifs = known_motifs + tmp
	print 'xml parse time', time.ctime()

    new_motifs = None
    # handle new motifs
    if opts.denovo:
        print 'Performing Denovo'
        # new motifs directory and file name
        #if opts.new_motifs_directory[-1] != '/':
        #    opts.new_motifs_directory += '/'
        if not os.path.exists(opts.new_motifs_directory):
            os.makedirs(opts.new_motifs_directory)
        new_motifs_file_name = os.path.join(opts.new_motifs_directory,
                                            opts.new_motifs)
        # find new motifs in ChIP regions
        new_motifs = chip_regions.mdmodule(width=width)
        
        # save new motifs in XML file
        new_motifs_file = open(new_motifs_file_name, 'w')
        new_motifs_file.write(new_motifs.to_xml())
        new_motifs_file.close()

    # set all_motifs to either: new_motifs, known_motifs or new + known motifs
    if new_motifs is None and known_motifs is None:
        #Error--but for now, just an empty thing
        all_motifs = None
    elif new_motifs is None:
        all_motifs = known_motifs
    elif known_motifs is None:
        all_motifs = new_motifs
    else:
        all_motifs = new_motifs + known_motifs

    # run seqpos statistic for each motif
    print "running seqpos stat for each motif..."
    all_motifs.seqpos(chip_regions, width=width)

    #cull the motifs
    sig_motifs = all_motifs.cull(pvalcutoff, maxmotif)

    print len(sig_motifs),"out of",len(all_motifs),"are significantly near peak centers."

    #filter by species?
    if opts.species_list is not None:
        species_list = opts.species_list.split(',')
        sig_motifs = sig_motifs.filterBySpecies(species_list)
    
    if sig_motifs is not None:
        # set motif to cluster on ( source motif from database or observed motif )
        sig_motifs.set_cluster_motif(motif_type='observed')
        # cluster all motifs using hierarchical clustering
        if len(sig_motifs) > 0: 
            motif_tree = sig_motifs.cluster()
            print 'clustering time', time.ctime()
        else: #No motifs found! return an empty MotifTree
            print 'No motifs found'
            motif_tree = MotifTree(None, None, None)
        # using SeqPos, score all motifs in motif tree against the ChIP regions
        #motif_tree.seqpos(chip_regions)
	#print 'seqpos time', time.ctime()

        # MOTIF Analysis complete, prepare output
        # output directory and file name
        #if opts.output_directory[-1] != '/':
        #    opts.output_directory += '/'
        if not os.path.exists(opts.output_directory):
            os.makedirs(opts.output_directory)
        output_file_name = os.path.join(opts.output_directory, opts.output)
        # images directory
        img_dir = opts.img_output_directory
        if not os.path.exists(img_dir):
            os.makedirs(img_dir)
        
        # save motif clustering to HTML file
        output_file = open(output_file_name, 'w')
        output_file.write(motif_tree.to_html(dst_dir=opts.output_directory,
                                             img_dir=img_dir))
        output_file.close()
	print 'html time', time.ctime()

        # copy static images to images directory
        for filename in os.listdir(STATIC_IMG_DIR):
            if os.path.isfile(filename):
                shutil.copy(STATIC_IMG_DIR + filename, img_dir)
