#! /usr/bin/env python

"""A module for sequence motifs."""

import os
import sys
import time
import math
import copy
import distutils.dir_util

import mdseqpos.settings as settings
import mdseqpos.count

import re

os.environ['DJANGO_SETTINGS_MODULE'] = 'mdseqpos.settings'

import numpy
import motiflogo

import xml.etree.ElementTree as etree
import xml.dom.minidom as minidom

from binarytree import BinaryTree
from django.template import loader, Context

# NOTE: SCHEMA_FILE_NAME is not used, I think I'll remove it soon
#SCHEMA_FILE_NAME = 'motif.xsd'

MUALPHA0 =  0.296 
MUALPHA1 =  0.641
SDALPHA0 =  0.943
SDALPHA1 = -0.111 

MOTIFMIN = 1e-3

import bayesian_motif_comp
from mdseqpos._seq import seqscan
import taolib.CoreLib.BasicStat.Prob as Prob

SENSE = bayesian_motif_comp.SENSE
ANTISENSE = bayesian_motif_comp.ANTISENSE

def revcomp(s):
    t = s[:]
    t.reverse()
    rc = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N' }
    t = [ rc[x] for x in t ]
    return t

class Motif:
    """A class for a sequence motif.
    
    Attributes:
    
    self.id : string, required
        A string containing a unique ID identifying the motif.  
    
    self.status : string, optional
        A string describing the status of the motif, such as "putative", 
        "confirmed", or "obsolete".
    
    self.source : string, optional
        A string describing the source of the motif, such as "Transfac" or 
        "PBM".
    
    self.sourcefile : string, optional
        A string containing the path to the source file of the motif.
    
    self.species : list of strings, optional
        A list of species names, such as ["mouse", "human", "yeast"].
    
    self.entrez : list of integers, optional
        A list of integers representing the Entrez Gene IDs of the
        transcription factor that bind to the motif, such as [6908, 2099].

    self.refseqs : list of strings, optional
        A list of the refseq strings, e.g. 'NM_001141919'

    self.symbols : list of strings, optionally empty
        A list of the transcription factors that bind to the motif, e.g.
        'ERa'.
    
    self.synonyms : list of strings, optional
        A list of synonyms for the transcription factor that binds to the
        motif, such as ['ESR', 'ESR1'].
    
    self.fullname : string, optional
        The full name of the transcription factor that binds to the motif,
        such as "estrogen receptor alpha".
    
    self.dbd : string, optional
        The full name of the DNA binding domain of the transcription factor
        that binds to the motif, such as "estrogen receptor region C".
    
    self.pmid : integer, optional
        An integer representing the Pubmed ID of the article in which the
        motif was published.
    
    self.pssm : NumPy float array (N rows, 4 columns), required
        A NumPy array containing the position specific scoring matrix (PSSM) 
        for the motif.  The array should contain N rows, each corresponding to 
        a single position in the PSSM, and 4 columns, each corresponding to
        the 4 DNA basepairs A, C, G, T, in order.  Each element in the matrix 
        should be a float, and each row of the matrix should sum to 1.0.
    
    self.numseqs : integer, optional
        The number of sequences used to derive the motif's position specific 
        scoring matrix (PSSM).
    
    self.curators : list of strings, optional
        A list of the full names of the curators, such as ["Cliff Meyer", 
        "Xiaole Shirley Liu", "Mayako Michino"].
    
    self.seqpos_results : dictionary
        A dictionary containing information about the SeqPos scores for the 
        motif.  The dictionary should contain the following keys:
          ### ADD THIS LATER ###
          ### (Miyako, the self.seqpos_results attribute isn't saved into the
               XML file, so you won't have to deal withthis. ###"""
    
    ##################
    # initialization #
    ##################
    
    def __init__(self, motif_id=None, motif_pssm=None):
        """Initialize the motif with PSSM m."""
        self.id = motif_id
        self.status = None
        self.source = None
        self.sourcefile = None
        self.species = []
        self.entrez = []
        self.refseqs = []
        self.symbols = []
        self.synonyms = []
        self.fullname = None
        self.dbd = None
        self.pmid = None
        if motif_pssm is None:
            self.pssm = None
        else:
            self.pssm = numpy.array(motif_pssm, float)
            if self.pssm.shape[1] != 4:
                raise ValueError, "motif PSSM must have 4 columns"
        self.numseqs = None
        self.curators = []
        self.seqpos_results = None
        # reverse complements need to be taken to make tree consistent
        self.antisense = False 
    
    def clean_str(self, attribute):
        """Clean up a string attribute."""
        value = getattr(self, attribute)
        if value is None:
            return
        string = str(value)
        if string.isspace():
            setattr(self, attribute, None)
        else:
            setattr(self, attribute, string)
    
    def clean_int(self, attribute):
        """Clean up an integer attribute."""
        value = getattr(self, attribute)
        if value is None:
            return
        try:
            integer = int(value)
            setattr(self, attribute, integer)
        except ValueError:
            print >> sys.stderr, "warning: cannot convert %s value '%s' to an integer" % (attribute, value)
            setattr(self, attribute, None)
    
    def clean_str_list(self, attribute):
        """Clean up a string list attribute."""
        value = getattr(self, attribute)
        if value is None:
            setattr(self, attribute, [])
            return
        if type(value) not in (list, tuple):
            value = [value]
        string_list = []
        for elem in value:
            if elem is None:
                continue
            string = str(elem)
            if not string.isspace():
                string_list.append(string)
        setattr(self, attribute, string_list)
    
    def clean_int_list(self, attribute):
        """Clean up an integer list attribute."""
        value = getattr(self, attribute)
        if attribute is None:
            setattr(self, attribute, [])
            return
        if type(value) not in (list, tuple):
            value = [value]
        integer_list = []
        for elem in value:
            if elem is None:
                continue
            try:
                integer = int(elem)
            except ValueError:
                print >> sys.stderr, "warning: cannot convert %s element '%s' to an integer" % (attribute, elem)
                continue
            integer_list.append(integer)
        setattr(self, attribute, integer_list)
    
    def cleanup(self):
        """Clean up formatting for all motif attributes."""
        for attribute in ('id', 'status', 'source', 'sourcefile', 'fullname', 'dbd'):
            self.clean_str(attribute)
        for attribute in ('pmid', 'numseqs'):
            self.clean_int(attribute)
        for attribute in ('species', 'refseqs', 'symbols', 'synonyms', 'curators'):
            self.clean_str_list(attribute)
        for attribute in ('entrez',):
            self.clean_int_list(attribute)
        if self.pssm is not None:
            try:
                self.pssm = numpy.array(self.pssm, float)
            except ValueError:
                print >> sys.stderr, "warning: cannot convert motif PSSM to a NumPy array:"
                print >> sys.stderr, self.pssm
                self.pssm = None
            if self.pssm.shape[1] != 4:
                print >> sys.stderr, "warning: motif PSSM must have 4 columns:"
                self.pssm = None
        
    ##################
    # static methods #
    ##################
 
    @staticmethod
    def isduplicate(motif1, motif2):
        """Return True if two motifs are the same, and False otherwise."""
        return motif1.id == motif2.id
 
    @staticmethod
    def merge(motif1, motif2, offset=0, antisense=False, merged_id=None):
        """Return a motif generated by merging two motifs."""
        merged = Motif()
        merged.id = merged_id
        merged.pssm = bayesian_motif_comp.merge(motif1.pssm, motif2.pssm, offset, antisense)
        merged.antisense = antisense
        return merged

    @staticmethod
    def merge_observed(motif1, motif2, offset=0, antisense=False, merged_id=None):
        """Return a motif generated by merging two motifs."""
        merged = Motif()
        merged.id = merged_id
        #merged.seqpos_results = { 'pssm': bayesian_motif_comp.merge( motif1.seqpos_results['pssm'], motif2.seqpos_results['pssm'], offset, antisense ) }
        merged.pssm_for_tree = bayesian_motif_comp.merge( motif1.pssm_for_tree, motif2.pssm_for_tree, offset, antisense )
        merged.antisense = antisense
        return merged
 
    @staticmethod
    def similarity(motif1, motif2, metric='Bayesian'):
        """Return a score for the similarity between two motifs.
        offset -- number of basepairs to shift the first motif
        antisense -- whether to take the reverse complement of the first motif"""
        similarity_score, offset, antisense = bayesian_motif_comp.BLiC_score(motif1.pssm, motif2.pssm)
        antisense = bool(antisense)
        return similarity_score, offset, antisense

    @staticmethod
    def similarity_observed(motif1, motif2, metric='Bayesian'):
        """Return a score for the similarity between two motifs.
        
        offset -- number of basepairs to shift the first motif
        antisense -- whether to take the reverse complement of the first motif"""
        #similarity_score, offset, antisense = bayesian_motif_comp.BLiC_score( motif1.seqpos_results['pssm'], motif2.seqpos_results['pssm'] )
        similarity_score, offset, antisense = bayesian_motif_comp.BLiC_score( motif1.pssm_for_tree, motif2.pssm_for_tree )
        antisense = bool(antisense)
        return similarity_score, offset, antisense

 
    #################
    # input methods #
    #################
    def pre_process_file(self, input_txt):
        """We need to pre-process the pssm input to ensure that it is in
        format described in 'from_flat_file'.  This may mean converting
        certain string literals to characters.  For example, galaxy converts
        [ --> __ob__ and ] --> to __cb__.  This method will convert them
        back"""
        d = [("__ob__", "["),
                ("__cb__", "]")]
        for s in d:
            input_txt = input_txt.replace(s[0], s[1])
        return input_txt
    
    def from_flat_file(self, filename):
        """Read in a motif from a text file which describes the pssm- the
        pssm is in the following format, example:
        [[0.1, 0.1, 0.7, 0.1],
         [0.1, 0.7, 0.1, 0.1],
         [0.7, 0.1, 0.1, 0.1],
         [0.1, 0.1, 0.1, 0.7]]
         """
        f = open(filename)
        pssm_txt = "\n".join([line for line in f]) #flatten file into a string
        pssm_txt = self.pre_process_file(pssm_txt) #translate weird literals
        
        # use regular expression to pick out the rows
        row_pattern = '\[\s*(0\.\d+)\s*\,\s*(0\.\d+)\s*\,\s*(0\.\d+)\s*\,\s*(0\.\d+)\s*\]'
        #function that converts a tuple of float strings to float objs
        convert_tpl = lambda x: [float(i) for i in x]
        pssm = [convert_tpl(r) for r in re.findall(row_pattern, pssm_txt)]

        #The pssm is a numpy.ndarray type, so we need to convert it
        self.pssm = numpy.array(pssm)
    
    def from_transfac(self):
        """Import motif from Transfac format."""
        pass
    
    def from_xml_file(self, xml_file_name):
        """Read list of motifs from XML file."""
        # parse motif XML
        tree = etree.parse(xml_file_name)
        # retrieve root tag
        root = tree.getroot()
        if root.tag != 'motif':
            raise ValueError, "root element of XML should be <motif>"
        # extract motif data
        self.from_xml_etree(root)

    def correct_pssm(self):
        """NOTE: this method is a hack.  In the databases, some pssms
        contain probability values that are absolutes, e.g. 0.0 or 1.0.
        We need to adjust 0 values to be 0.01 (b/c we stop at just 100
        motifs); and 1 -> 0.97.
        For rows where there are just there's just one or two 0s, we
        need to decrement the other values in the row accordingly.

        Additional note: i was told by tao, that we need only check for
        the case where there is 1.0 and three 0s.  So that would make
        this algorithm obsolete.
        """

        for row in self.pssm:
            numZeros = 0 #first pass, count the number of 0s
            for val in row:
                if val == 0:
                    numZeros += 1
            if numZeros > 0: #need to modify this row
                #decr is the amount we reduce the NON-zero values by
                #e.g. if there is one 0.0, then we need to decr
                #by 0.01/3 = 0.003333
                decr = numZeros*0.01/(4 - numZeros)
                ##print "OLD: %s" % row
                #make the corrections for the row
                for i, prob in enumerate(row):
                    if prob == 0: 
                        row[i] = 0.01
                    else:
                        row[i] -= decr
                ##print "New: %s %s %s" % (row, decr, numZeros)
    
    def from_xml_etree(self, root):
        """Read motif from an XML ElementTree.
        
        Any data that was stored in Motif object will be overwritten."""
        assert root.tag == 'motif'
        # retrieve motif ID
        self.id = root.attrib['id']
        # retrieve status
        self.status = root.find('status')
        if self.status is not None:
            self.status = self.status.text
        # retrieve source
        self.source = root.find('source')
        if self.source is not None:
            self.source = self.source.text
        # retrieve path to source file
        self.sourcefile = root.find('sourcefile')
        if self.sourcefile is not None:
            self.sourcefile = self.sourcefile.text
        # retrieve list of species
        specieslist = root.find('specieslist')
        if specieslist is not None:
            self.species = [species.text.strip() for species in specieslist.findall('species')]
        # retrieve Entrez GeneID
        self.entrezlist = root.find('entrezlist')
        if self.entrezlist is not None:
            self.entrez = [int(entrez.text) for entrez in self.entrezlist.findall('entrez')]
        # retrieve refseq list
        self.refseqlist = root.find('refseqlist')
        if self.refseqlist is not None:
            refs = self.refseqlist.findall('refseq')
            self.refseqs = [r.text.strip() for r in refs]
        
        # retrieve list of symbols, ie. factors
        # NOTE: its weird that self.symbols is an XML elm, and then a list
        # but i'm just going to follow the format
        self.symbols = root.find('symbollist')
        if self.symbols is not None:
            syms = self.symbols.findall('symbol')
            self.symbols = [factor.text.strip() for factor in syms]
        
        # retrieve list of synonyms
        self.synonyms = [synonym.text for synonym in root.findall('synonym')]
        # retrieve full name
        self.fullname = root.find('fullname')
        if self.fullname is not None:
            self.fullname = self.fullname.text
        # retrieve DNA binding domain
        self.dbd = root.find('dbd')
        if self.dbd is not None:
            self.dbd = self.dbd.text
        # retrieve PubMed ID
        self.pmid = root.find('pmid')
        if self.pmid is not None:
            self.pmid = int(self.pmid.text)
        # retrieve PSSM
        self.pssm = root.find('pssm')
        if self.pssm is not None:
            pos_elements = self.pssm.findall('pos')
            self.pssm = numpy.empty((len(pos_elements), 4), float)
            for i, pos in enumerate(pos_elements):
                assert int(pos.attrib['num']) == i + 1
                self.pssm[i,:] = [float(base.text) for base in pos]

        # Correct the PSSM values - see the method descr. above - OFF for now
        #self.correct_pssm()
        
        # retrieve number of sequences
        self.numseqs = root.find('numseqs')
        if self.numseqs is not None:
            self.numseqs = int(self.numseqs.text)
        # retrieve list of curators
        self.curators = [curator.text for curator in root.findall('curator')]
        # clear seqpos results
        self.seqpos_results = None
        #print self
    
    ##################
    # output methods #
    ##################
    
    def __str__(self):
        """Generate string representation of motif."""
        motif_string = []
        motif_string.append("Motif ID: %s" % self.id)
        motif_string.append("Status: %s" % self.status)
        motif_string.append("Source: %s" % self.source)
        motif_string.append("Source File: %s" % self.sourcefile)
        if len(self.species) == 0:
            motif_string.append("Species: None")
        else:
            motif_string.append("Species: %s" % ', '.join(self.species))
        if len(self.entrez) == 0:
            motif_string.append("Entrez GeneID: None")
        else:
            motif_string.append("Entrez GeneID: %s" % ', '.join(str(entrez) for entrez in self.entrez))
            
        if len(self.refseqs) == 0:
            motif_string.append("Refseqs: None")
        else:
            motif_string.append("Refseqs: %s" % ', '.join(str(ref) for ref in self.refseqs))
            
        if self.symbols is None or len(self.symbols) == 0:
            motif_string.append('Symbols: None')
        else:
            motif_string.append("Symbols: %s" % ','.join(self.symbols))
            
        if len(self.synonyms) == 0:
            motif_string.append("Synonyms: None")
        else:
            motif_string.append("Synonyms: %s" % ', '.join(self.synonyms))
        motif_string.append("Full Name: %s" % self.fullname)
        motif_string.append("DNA Binding Domain: %s" % self.dbd)
        motif_string.append("PMID: %s" % self.pmid)
        if self.pssm is None:
            motif_string.append("PSSM: None")
        else:
            motif_string.append("PSSM:        A      C      G      T  ")
            for i, pos in enumerate(self.pssm):
                motif_string.append(" |%6d   %4.2f   %4.2f   %4.2f   %4.2f" % (i + 1, pos[0], pos[1], pos[2], pos[3]))
        motif_string.append("Number of Sequences: %s" % self.numseqs)
        if len(self.curators) == 0:
            motif_string.append("Curators: None")
        else:
            motif_string.append("Curators: %s" % ', '.join(self.curators))
        if self.seqpos_results is None:
            motif_string.append("SeqPos Results: None")
        else:
            motif_string.append("SeqPos Results:")
            motif_string.append(" |  Number of Hits ... %d" % self.seqpos_results['numhits'])
            motif_string.append(" |  Cutoff ........... %.2f" % self.seqpos_results['cutoff'])
            motif_string.append(" |  Z-score .......... %.2f" % self.seqpos_results['zscore'])
            motif_string.append(" |  p-value .......... %.2f" % self.seqpos_results['pvalue'])
        motif_string = '\n'.join(motif_string)
        return motif_string

    def to_xml(self, pretty_print=False):
        """Return motif as XML."""
        root = self.to_xml_etree()
        if root is not None:
            motif_as_xml = etree.tostring(root)
            if pretty_print is True:
                motif_as_xml = minidom.parseString(motif_as_xml).toprettyxml(indent="  ")
            return motif_as_xml
    
    def to_xml_etree(self):
        """Return motif as an XML ElementTree."""
        # motif ID and motif PSSM must be defined
        if self.id is None:
            print >> sys.stderr, "motif ID 'self.id' must be defined for conversion to XML"
            return
        if self.pssm is None:
            print >> sys.stderr, "motif PSSM 'self.pssm' must be defined for conversion to XML"
            return
        # root element
        root = etree.Element('motif')
        root.attrib['id'] = str(self.id)
        # status element
        if self.status is not None:
            status_element = etree.SubElement(root, 'status')
            status_element.text = str(self.status)
        # source element
        if self.source is not None:
            source_element = etree.SubElement(root, 'source')
            source_element.text = str(self.source)
        # sourcefile element
        if self.sourcefile is not None:
            sourcefile_element = etree.SubElement(root, 'sourcefile')
            sourcefile_element.text = str(self.sourcefile)
        # species elements
        for species in self.species:
            species_element = etree.SubElement(root, 'species')
            species_element.text = str(species)
        # entrez element (parent: entrezlist, children: entrez)
        entrezlist = etree.SubElement(root, 'entrezlist')
        for entrez in self.entrez:
            entrez_element = etree.SubElement(entrezlist, 'entrez')
            entrez_element.text = str(entrez)
        # refseqs element (parent: refseqlist, children: refseq)
        refseqlist = etree.SubElement(root, 'refseqlist')
        for refseq in self.refseqs:
            refseq_element = etree.SubElement(refseqlist, 'refseq')
            refseq_element.text = refseq
        # symbols element- (parent: symbollist, children: symbol)
        symbollist = etree.SubElement(root, 'symbollist')
        for sym in self.symbols:
            sym_elm = etree.SubElement(symbollist, 'symbol')
            sym_elm.text = sym
        # synonym elements
        for synonym in self.synonyms:
            synonym_element = etree.SubElement(root, 'synonym')
            synonym_element.text = str(synonym)
        # fullname element
        if self.fullname is not None:
            fullname_element = etree.SubElement(root, 'fullname')
            fullname_element.text = str(self.fullname)
        # dbd element
        if self.dbd is not None:
            dbd_element = etree.SubElement(root, 'dbd')
            dbd_element.text = str(self.dbd)
        # pmid element
        if self.pmid is not None:
            pmid_element = etree.SubElement(root, 'pmid')
            pmid_element.text = str(self.pmid)
        # pssm element
        pssm_element = etree.SubElement(root, 'pssm')
        for i, pos in enumerate(self.pssm):
            pos_element = etree.SubElement(pssm_element, 'pos')
            pos_element.attrib['num'] = str(i + 1)
            A_element = etree.SubElement(pos_element, 'A')
            A_element.text = "%.2f" % pos[0]
            C_element = etree.SubElement(pos_element, 'C')
            C_element.text = "%.2f" % pos[1]
            G_element = etree.SubElement(pos_element, 'G')
            G_element.text = "%.2f" % pos[2]
            T_element = etree.SubElement(pos_element, 'T')
            T_element.text = "%.2f" % pos[3]
        # numseqs element
        if self.numseqs is not None:
            numseqs_element = etree.SubElement(root, 'numseqs')
            numseqs_element.text = str(self.numseqs)
        # curator elements
        for curator in self.curators:
            curator_element = etree.SubElement(root, 'curator')
            curator_element.text = str(curator)
        return root
    
    def logo(self, width, height, logo_dir="img/"):
        """Generate motif logo and return its file name."""
        # create directory for holding motif logos if directory does not
        # already exist
        if not os.path.exists(logo_dir):
            os.makedirs(logo_dir)
        # generate temporary FASTA file with 1000 sequences, with base
        # frequencies approximating the motif PSSM
        fasta_file_name = '%stemp.fa' % logo_dir
        motiflogo.fasta(self.pssm_for_tree, fasta_file_name, 1000)
        # generate motif logo using the temporary FASTA file
        #logo_file_name = '%s%s_%dx%d' % (logo_dir, self.id, width, height)
        logo_name = '%s_%dx%d' % (self.id, width, height)
        logo_file_name = os.path.join(logo_dir, logo_name)
        command = "%s -f %s -o %s -w %d -h %d -F PNG -c -n -Y" % (os.path.join(settings.DEPLOY_DIR, 'weblogo', 'seqlogo'), fasta_file_name, logo_file_name, width / 20.0, height / 20.0)
        os.system(command)
        # delete the temporary FASTA file
        os.remove(fasta_file_name)
        # return the filename of the motif logo
        logo_file_name += '.png'
        return logo_file_name
 
    
    def to_html(self, img_dir="img/"):
        """Return motif as pretty HTML for display purposes."""
        # generate a motif logo
        logo_file_name = self.logo(192, 120, logo_dir=img_dir)
        # prepare the motif PSSM for display
        motif_pssm = []
        if self.pssm is None:
            motif_pssm.append("PSSM: None")
        else:
            motif_pssm.append("        A      C      G      T  ")
            for i, pos in enumerate(self.pssm_for_tree):
                motif_pssm.append("%3d   %4.2f   %4.2f   %4.2f   %4.2f" % (i + 1, pos[0], pos[1], pos[2], pos[3]))
        # prepare SeqPos results for display
        seqpos_results = []
        if self.seqpos_results is None or ( 'numhits' not in self.seqpos_results ):
            seqpos_results.append("SeqPos Results: None")
        else:
            seqpos_results.append("SeqPos Results:")
            seqpos_results.append(" |  Number of Hits ... %d" % self.seqpos_results['numhits'])
            seqpos_results.append(" |  Cutoff ........... %.2f" % self.seqpos_results['cutoff'])
            seqpos_results.append(" |  Z-score .......... %.2f" % self.seqpos_results['zscore'])
            seqpos_results.append(" |  p-value .......... %.4e" % self.seqpos_results['pvalue'])
        # load Django template for motifs
        t = loader.get_template('motif.html')
        # load context for Django template
        c = Context({
            'img_dir': img_dir,
            'motif_name': self.id,
            'motif_logo': logo_file_name,
            'motif_pssm': '\n'.join(motif_pssm),
            'seqpos_results': '\n'.join(seqpos_results),
        })
        motif_as_html = t.render(c)
        return motif_as_html

    # 08-23-10 NOTE: I think this pssm_to_json function can be a lot prettier!!
    def pssm_to_json(self):
        """Returns the motif matrix as a valid JSON matrix"""
        if self.pssm_for_tree is None:
            return "[[]]" #empty matrix

        tmp = "[" #open mat
        for i,row in enumerate(self.pssm_for_tree):
            tmp += "[" #open row
            for j, val in enumerate(row):
                if j != (len(row) - 1): #handle the commas
                    tmp += ("%.4f" % val)+", " #to 4 significant digits
                else:
                    tmp += ("%.4f" % val)
            if i != (len(self.pssm_for_tree) -1):
                #handle commas while closing row
                tmp += "],"
            else:
                tmp += "]"
        tmp +="]" #close mat
        return tmp            
        
    def to_json(self, dst_dir="results/", img_dir="img/"):
        """Returns the json representation of the motif"""
        #NOTE: the A if X else B,isn't supported in python2.4
        if self.id is not None:
            id = "'id':'"+self.id+"'"
        else:
            id = "'id':'None'"            
        if self.symbols is not None and len(self.symbols) > 0:
            factors = "'factors':[%s]" % ','.join(["'%s'" % s for s in self.symbols])
        else:
            factors = "'factors':[]"

        if self.entrez is not None and len(self.entrez) > 0:
            entrez = "'entrezs':[%s]" % ','.join(["'%s'" % s for s in self.entrez])
        else:
            entrez = "'entrezs':[]"

        if self.refseqs is not None and len(self.refseqs) > 0:
            refseqs = "'refseqs':[%s]" % ','.join(["'%s'" % s for s in self.refseqs])
        else:
            refseqs = "'refseqs':[]"

        if self.species is not None and len(self.species) > 0:
            species = "'species':%s" % map(lambda s:'%s' % s.strip(),
                                           self.species)
        else:
            species = "'species':[]"
       
        if self.seqpos_results == None:
            hits = "'None'"
            cutoff = "'None'"
            zscore = "'None'"
            pvalue = "'None'"
            mu = "'None'"
        else:
            hits   = str(self.seqpos_results['numhits'])
            cutoff = "%.4f" % self.seqpos_results['cutoff']
            zscore = "%.4f" % self.seqpos_results['zscore']
            mu     = "%.4f" % self.seqpos_results['meanposition']
            pvalue = "%.4e" % self.seqpos_results['pvalue']


        logoImg = self.logo(192, 120, logo_dir=img_dir)
        #since it is a relative path, parse out the dst_dir, from logoImg path
        logoImg = logoImg.replace(dst_dir, "")
        
        return "{"+id+", "+factors+','+entrez+','+refseqs+','+species+\
               ", 'consensus':'None', 'pssm':"+self.pssm_to_json()+\
               ", 'logoImg':'"+logoImg+"', 'hits':"+hits+", 'cutoff':"+cutoff+\
               ", 'zscore':"+zscore+", 'pval':"+pvalue+", 'position':"+mu+"}"
    
    ####################
    # analysis methods #
    ####################
    def seqpos_stat(self, length, start, end, score, offset):
        """BINOCH: FROM binoch/src/motifanalysis/nsd.py, a method to calculate
        motif stats"""
        
        start   = numpy.array(start,float)
        end     = numpy.array(end,float)
        length  = numpy.array(length,float)
        #print "Start is:\n%s\n" % start
        #print "End is:\n%s\n" % end
        #print "Length is:\n%s\n" % length
        #print "Score is:\n%s\n" % score
        #print "Offset is:\n%s\n" % offset
        
        MINHITS = 20
        score   = -score
        frac    = numpy.abs( ( start + end - offset )/(length-offset) - 1.0 )
        #print "frac is:\n%s\n" % frac

        #NOTE: argsort returns the list of indices sorted by score values,
        #e.g. idx[0] is the index of the largest score (b/c we above
        #score=-score, and score[idx[0]] is largest score
        #try this instead:
        #idx = [i for (i,s) in sorted(enumerate(score), lambda (i,s): s, reverse)]
        idx     = score.argsort()
        mu_bias = 0.5
        #print "idx is:\n%s\n" % idx
        
        cumfrac, meanfracpos_t, zscore_t = 0,[],[]
        zscore_min, kmin, numhits, cuthits, meanpos = 99, -1, 0, 0, 0

        #I BELIEVE this is a bug: numhits should not be some enumerated index!
        for numhits,i in enumerate(idx):
            #print "numhits,i is (%s,%s)" % (numhits, i)
            
            if frac[i] > 1.0 or frac[i] < 0.0:
                #print "frac out of bounds",start[i],end[i], length[i], frac[i]
                raise(ValueError,"Fraction out of bounds %4.2f" % frac[i] )

            cumfrac += frac[i]
            meanfracpos_t.append( cumfrac/numhits )
            if numhits > MINHITS:
                zscore = ( cumfrac - numhits*mu_bias )/( math.sqrt( numhits/12.0 ) )
                #print 'zzz', numhits,score[i],frac[i],meanfracpos_t[numhits],zscore
                if zscore < zscore_min:
                    zscore_min = zscore
                    kmin = -score[i]
                    cuthits = numhits
                    meanpos = cumfrac/numhits - 0.5

        if cuthits > MINHITS:
            # NOTE: w/ the pvalue, we are assuming the normal distribution
            # get the raw p-value score
            z_mean_adj = - (MUALPHA0 + MUALPHA1*math.log(math.log(numhits)))
            z_sd_adj = SDALPHA0 + SDALPHA1*math.log(math.log(numhits))
            zscore_min_adjusted = ( zscore_min - z_mean_adj )/z_sd_adj  
            pvalue = max(1e-30, Prob.normal_cdf(zscore_min_adjusted))
            #print "#----->", cuthits, z_mean_adj, z_sd_adj, zscore_min_adjusted, pvalue
        else:
            zscore_min_adjusted = 0
            pvalue = 1.0
            #M = 0.25*numpy.ones( self.pssm.shape )
            
        results = {
            'numhits': cuthits,
            'meanposition': meanpos,
            'cutoff': kmin,
            'zscore': zscore_min_adjusted,
            'pvalue': pvalue,
            #'pssm': M
            }
        #print results
        return results
    
    def seqpos(self, chip_regions, preprocessed_regions=False, width=600, margin=50 ):
        """Score motif on how centrally located they are within each ChIP region.
        
        ChIP regions should be given as a ChipRegions object.
        The results of SeqPos will be stored as properties of self.
 
        Adapted from Cliff Meyer's ChIP_region.CentralMotifScan() method."""

        # run SeqPos magic
        # seqscan(regions.sequence, self.pssm)
        # ...
        #LEN NOTE: width is a magic variable--it should be defined somewhere,
        #but it isn't, and i'm just using my best guess as to what it should be
        #set to (after consulting with Cliff)
	# pad sequences with margin so that hits can be detected at the extremities 

        # prepare sequences
        if preprocessed_regions is False:
            regions = chip_regions.copy()
            regions.unique()
            regions.trim(width + margin )
            regions.read_sequence()
            #regions.ext_read_sequence()
        else:
            regions = chip_regions

        # LEN: BINOCH UPGRADE
        bgseqprob_mat = mdseqpos.count.count(regions.sequence)
        markov_order = 2
        prob_option = mdseqpos._seq.MAX_OPTION

        #print "Regions sequence: \n%s\n" % regions.sequence
        
        s_idx, start, end, orient, score = seqscan(regions.sequence, self.pssm,
                                                   bgseqprob_mat, markov_order,
                                                   prob_option)
        lengths = [len(regions.sequence[s_idx[elm]]) for elm in s_idx]

        #adjust score
        score = (numpy.log(score + MOTIFMIN))
        self.seqpos_results = self.seqpos_stat(numpy.array(lengths), start,
                                               end, score, len(self.pssm))
        
        # retrieve sequences above cutoff #start, end
        #fracpos = (abs(0.5*(locs_t[:,1] + locs_t[:,2]) - \
        #               (margin+width)/2.0)) / (width/2.0)
        fracpos = (abs(0.5*(start + end) - (margin + width)/2.0)) / (width/2.0)
        #print "FRACPOS is: %s" % fracpos
        seq = []
        for j,elem2 in enumerate(fracpos):
            if elem2 <= 1.0:
                t = list(regions.sequence[int(s_idx[j])])
                t = t[int(start[j]):int(end[j])]
                if orient[j] == ANTISENSE:
                    seq.append(revcomp(t))
                else:
                    seq.append(t)

        #print "seq is: %s" % seq
        T = numpy.array(seq)
        M = []
        if (len(T) > 0): #if its a non-empty sequence!
            for elem in ['A','C','G','T']:
                M.append( list( (1*( T == elem )).sum(axis=0) ) )
            #print "M1 is %s" % M
            M = numpy.array( M, float )
            M = M.transpose()
            #print "M2 is %s" % M
            for i in range(M.shape[0]):
                M[i,:] = M[i,:]/(M[i,:].sum())
            #print "M3 is %s" % M
        #print "PSSM is: %s" % self.pssm
        self.seqpos_results['pssm'] = M
        # LEN: END BINOCH UPGRADE        

class MotifList(list):
    """A class for a list of sequence motifs."""
    
    ##################
    # initialization #
    ##################
    
    def __init__(self, iterable=[]):
        list.__init__(self, iterable)
        self.check_type()
    
    ######################################################
    # inherited methods, supplemented with type checking #
    ######################################################
    
    def __add__(self, another_motif_list):
        new_motif_list = MotifList()
        new_motif_list.extend(self)
        new_motif_list.extend(another_motif_list)
        return new_motif_list
    
    def __setitem__(self, index, value):
        list.__setitem__(self, index, value)
        self.check_type()
    
    def __setslice__(self, index_start, index_end, iterable):
        list.__setslice__(self, index_start, index_end, iterable)
        self.check_type()
    
    def append(self, value):
        list.append(self, value)
        self.check_type()
    
    def extend(self, iterable):
        list.extend(self, iterable)
        self.check_type()
    
    def insert(self, index, value):
        list.append(self, index, value)
        self.check_type()

    def cull(self, pvalcutoff=0.001, maxmotifs=100):
        """Filter the list based on the pvalue and the max number
        of motifs allowed in the list"""
        sig_motifs = filter(lambda m: m.seqpos_results['pvalue'] <= pvalcutoff,
                            self)
        if len(sig_motifs) > maxmotifs:
            #sort and return top maxmotifs
            sorted_motifs = sorted(self, key = lambda elem: elem.seqpos_results['pvalue'])
            del sorted_motifs[maxmotifs:]
            return MotifList(sorted_motifs)
        else:
            return MotifList(sig_motifs)

    def filterBySpecies(self, species_list):
        """Filter the list by the species in species_list"""
        mapSpecies = {'hs' : 'Homo sapiens', 'mm' : 'Mus musculus',
                      'dm' : 'Drosophila melanogaster',
                      'ce' : 'Caenorhabditis elegans'}
        #MAP short names to long names
        sl = map(lambda s: mapSpecies[s], species_list)
        motifs = []
        for m in self:
            #take the motif IFF species = [], OR an element in species is
            #in species_list
            if m.species:
                found = False
                i = 0
                while (i < len(m.species) and not found):
                    if (m.species[i] in sl): found = True
                    i += 1
                if found: motifs.append(m)
            else:
                motifs.append(m)
        return MotifList(motifs)
    
    ######################
    # consistency checks #
    ######################
    
    def check_type(self):
        for elem in self:
            if not isinstance(elem, Motif):
                raise TypeError, "element %s is not of %s" % (elem, Motif)
    
    def check_unique(self):
        for i in xrange(len(self)):
            for j in xrange(i + 1, len(self)):
                if Motif.isduplicate(self[i], self[j]):
                    raise ValueError, "motifs %s and %s of %s are the same" % (i, j, self)
    
    #################
    # input methods #
    #################
    
    def from_xml_file(self, xml_file_name):
        """Read list of motifs from XML file and append to MotifList."""
        # parse motif XML
        tree = etree.parse(xml_file_name)
        # retrieve root tag
        root = tree.getroot()
        if root.tag != 'motifs':
            return  # error
        # extract motif data
        for child in root:
            new_motif = Motif()
            new_motif.from_xml_etree(child)
            self.append(new_motif)
    
    ##################
    # output methods #
    ##################
    
    def to_xml(self, pretty_print=False):
        """Write list of motifs to XML file."""
        # root <motifs> element
        root = etree.Element('motifs')
        # child <motif> elements
        for motif in self:
            motif_element = motif.to_xml_etree()
            root.append(motif_element)
        motifs_as_xml = etree.tostring(root)
        if pretty_print is True:
            motifs_as_xml = minidom.parseString(motifs_as_xml).toprettyxml(indent="  ")
        return motifs_as_xml
    
    ####################
    # analysis methods #
    ####################
    
    def seqpos(self, chip_regions, width=600, margin=50 ):
        """Score motifs on how centrally located they are within each ChIP
        region.
  
        ChIP regions should be given as a ChipRegions object.
        The results of SeqPos will be stored as properties each Motif object.
        
        Adapted from Cliff Meyer's ChIP_region.CentralMotifScan() method."""
        # preprocess regions by preparing sequences
        regions = chip_regions.copy()
        regions.unique()
        regions.trim(width+margin)
        #print "READ SEQUENCE START %s\n" % time.time()*1000
        regions.read_sequence()
        #regions.ext_read_sequence()
        #print "READ SEQUENCE END %s\n" % time.time()*1000
        # run SeqPos for each motif against the preprocessed regions
        for motif in self:
            motif.seqpos(regions, preprocessed_regions=True, width=width, margin=margin )

    def cluster(self):
        """Cluster list of motifs using hierchical clustering."""
	t = MotifTreeList(self).cluster()
        t.append_source_motifs()
        t.correct_motif_orientation( antisense = False )
        return t 

    def set_cluster_motif(self,motif_type='observed'):
        for elem in self:
            if motif_type == 'observed':
                #print elem.id
                elem.pssm_for_tree = copy.deepcopy(elem.seqpos_results['pssm'])
            if motif_type == 'source':
                elem.pssm_for_tree = copy.deepcopy(elem.pssm)

class MotifTree(BinaryTree):
    """A class for a tree of sequence motifs."""
 
    #######################
    #   initialization    #
    #######################
    
    def __init__(self, v=None, l=None, r=None):
        BinaryTree.__init__(self, v, l, r)
        self.check_type(Motif)
 
    ######################
    # motif manipulation # 
    ###################### 
    def append_source_motifs(self):
        """
        to the leaf nodes of the tree append the original pssm
        """
        if self.lchild:
            self.lchild.append_source_motifs() 
        if self.rchild:
            self.rchild.append_source_motifs() 
        if ( not self.lchild ) and ( not self.rchild ):
            self.lchild = copy.deepcopy(self)
            self.lchild.value.pssm_for_tree = copy.deepcopy(self.lchild.value.pssm)
            self.value.id = self.value.id + "_observed"
 
    def correct_motif_orientation( self, antisense = False):
        """
        antisense refers to whether a node was flipped on merging
        """

        #print self.value.id, antisense, self.value.antisense

        if antisense == True:

            try:
                self.value.pssm_for_tree = self.value.pssm_for_tree[::-1,::-1]
            except:
                print 'exception in motif reversal'            

            if self.value.antisense == True:
                if self.lchild:
                    self.lchild.correct_motif_orientation(antisense=False) 
                if self.rchild:
                    self.rchild.correct_motif_orientation(antisense=True)
            else:
                if self.lchild:
                    self.lchild.correct_motif_orientation(antisense=True) 
                if self.rchild:
                    self.rchild.correct_motif_orientation(antisense=True)

            self.value.antisense = False

        else:
  
            if self.value.antisense == False:
                if self.lchild:
                    self.lchild.correct_motif_orientation(antisense=False) 
                if self.rchild:
                    self.rchild.correct_motif_orientation(antisense=False)
            else:
                if self.lchild:
                    self.lchild.correct_motif_orientation(antisense=True) 
                if self.rchild:
                    self.rchild.correct_motif_orientation(antisense=False)
 
    ######################
    #   output methods   #
    ######################
 
    def to_xml(self):
        """Convert tree of motifs into a flat list of motifs in XML format."""
        return MotifList(self.traverse()).to_xml()
    
    def to_html(self, dst_dir='results/', img_dir='img/'):
        """Renders the MotifTree into HTML form, and copies all of the
        supporting files associated with that HTML output"""
        #copy over supporting files--i.e. everything in the template dir/static
        distutils.dir_util.copy_tree(os.path.join(settings.DEPLOY_DIR,
                                                  "django/", "static/"),
                                     dst_dir, update=1)
        return BinaryTree.to_new_html(self, "mdseqpos_out.html", dst_dir, img_dir)
    
    ####################
    # analysis methods #
    ####################
    
    def seqpos(self, chip_regions, width=600, margin=50 ):
        """Score motifs on how centrally located they are within each ChIP region.
        
        ChIP regions should be given as a ChipRegions object.
        The results of SeqPos will be stored as properties each Motif object.
        
        Adapted from Cliff Meyer's ChIP_region.CentralMotifScan() method."""
        # preprocess regions by preparing sequences
        regions = chip_regions.copy()
        regions.unique()
        regions.trim(width+margin)
        regions.read_sequence()
        #regions.ext_read_sequence()
        # run SeqPos for each motif against the preprocessed regions
        for motif in self.traverse():
            motif.seqpos(regions, preprocessed_regions=True, width=width, margin=margin)



class MotifTreeList(list):
    """A class for a list of motif trees."""
 
    ##################
    # initialization #
    ##################
 
    def __init__(self, iterable=[]):
        raw_list = []
        for elem in iterable:
            if isinstance(elem, MotifTree):
                raw_list.append(elem)
            elif isinstance(elem, Motif):
                raw_list.append(MotifTree(elem, None, None))
            else:
                raise TypeError, "element %s is not of %s" % (elem, MotifTree)
        list.__init__(self, raw_list)
        self.check_type()
        self.set_similarity()
 
    def set_similarity(self):
        # initialize similarity
        self.similarity = numpy.zeros((len(self), len(self)), float)
        self.offset = numpy.zeros((len(self), len(self)), int)
        self.antisense = numpy.zeros((len(self), len(self)), bool)
        # first, calculate upper triangle of the similarity matrix
        for i in xrange(len(self)):
            for j in xrange(i + 1, len(self)):
                print "computing similarity of %d and %d" % (i,j)
                self.similarity[i,j], self.offset[i,j], self.antisense[i,j] = Motif.similarity_observed(self[i].value, self[j].value)
        # then, use transpose of upper triangle for the lower triangle
        self.similarity += self.similarity.T
        self.offset += self.offset.T
        self.antisense += self.antisense.T
    
    ################################################################################
    # inherited methods, supplemented with type checking and similarity attributes #
    ################################################################################
    
    def __setitem__(self, index, value, prechecked_type=False):
        list.__setitem__(self, index, value)
        j = index
        for i in xrange(len(self)):
            if i != j:
                self.similarity[i,j], self.offset[i,j], self.antisense[i,j] = Motif.similarity(self[i].value, self[j].value)
                self.similarity[j,i] = self.similarity[i,j]
                self.offset[j,i] = self.offset[i,j]
                self.antisense[j,i] = self.antisense[i,j]
    
    def __setslice__(self, index_start, index_end, iterable):
        list.__setslice__(self, index_start, index_end, iterable)
        self.check_type()
        self.set_similarity()
    
    def append(self, value):
        list.append(self, value)
        self.check_type()
        self.similarity = numpy.insert(self.similarity, self.similarity.shape[0], 0.0, axis=0)
        self.similarity = numpy.insert(self.similarity, self.similarity.shape[1], 0.0, axis=1)
        self.offset = numpy.insert(self.offset, self.offset.shape[0], 0.0, axis=0)
        self.offset = numpy.insert(self.offset, self.offset.shape[1], 0.0, axis=1)
        self.antisense = numpy.insert(self.antisense, self.antisense.shape[0], 0.0, axis=0)
        self.antisense = numpy.insert(self.antisense, self.antisense.shape[1], 0.0, axis=1)
        j = len(self) - 1
        for i in xrange(len(self)):
            if i != j:
                self.similarity[i,j], self.offset[i,j], self.antisense[i,j] = Motif.similarity_observed(self[i].value, self[j].value)
                self.similarity[j,i] = self.similarity[i,j]
                self.offset[j,i] = self.offset[i,j]
                self.antisense[j,i] = self.antisense[i,j]
    
    def extend(self, iterable):
        list.extend(self, iterable)
        self.check_type()
        self.set_similarity()
    
    def insert(self, index, value):
        list.append(self, index, value)
        self.check_type()
        self.similarity = numpy.insert(self.similarity, index, 0.0, axis=0)
        self.similarity = numpy.insert(self.similarity, index, 0.0, axis=1)
        self.offset = numpy.insert(self.offset, self.offset.shape[0], 0.0, axis=0)
        self.offset = numpy.insert(self.offset, self.offset.shape[1], 0.0, axis=1)
        self.antisense = numpy.insert(self.antisense, self.antisense.shape[0], 0.0, axis=0)
        self.antisense = numpy.insert(self.antisense, self.antisense.shape[1], 0.0, axis=1)
        j = index
        for i in xrange(len(self)):
            if i != j:
                self.similarity[i,j], self.offset[i,j], self.antisense[i,j] = Motif.similarity(self[i].value, self[j].value)
                self.similarity[j,i] = self.similarity[i,j]
                self.offset[j,i] = self.offset[i,j]
                self.antisense[j,i] = self.antisense[i,j]
    
    def pop(self, index=None):
        if index is None:
            index = len(self) - 1
        list.pop(self, index)
        self.similarity = numpy.delete(self.similarity, index, axis=0)
        self.similarity = numpy.delete(self.similarity, index, axis=1)
        self.offset = numpy.delete(self.offset, index, axis=0)
        self.offset = numpy.delete(self.offset, index, axis=1)
        self.antisense = numpy.delete(self.antisense, index, axis=0)
        self.antisense = numpy.delete(self.antisense, index, axis=1)
    
    def remove(self, value):
        index = self.index(value)
        self.pop(index)
    
    def reverse(self):
        list.reverse(self)
        self.similarity = numpy.flipud(self.similarity)
        self.similarity = numpy.fliplr(self.similarity)
        self.offset = numpy.flipud(self.offset)
        self.offset = numpy.fliplr(self.offset)
        self.antisense = numpy.flipud(self.antisense)
        self.antisense = numpy.fliplr(self.antisense)
    
    def sort(self, cmp=None, key=None, reverse=False):
        list.sort(self, cmp, key, reverse)
        self.set_similarity()
    
    ######################
    # consistency checks #
    ######################
    
    def check_type(self):
        for i, elem in enumerate(self):
            if isinstance(elem, Motif):
                self[i] = MotifTree(elem, None, None)
            elif not isinstance(elem, MotifTree):
                raise TypeError, "element %s is not of %s" % (elem, MotifTree)
    
    ####################
    # analysis methods #
    ####################

    def cluster(self):
        infinity = 8000
        print "clustering... starting with %d trees in MotifTreeList" % len(self)
        while len(self) > 1:
            t1 = time.time()
            if len(self) > 2:
                #compute the similarity score of each motif compared to others
                #HEAVILY discount the diagonal of this matrix b/c we are only
                #interested in the score as compared to others, NOT to itself.
                #BUG: IF (k,k) is ever selected as the next merge, then
                #the binary tree collapses!!!
                similarity = self.similarity - (numpy.eye(len(self))*infinity)
                (i, j) = numpy.unravel_index(similarity.argmax(), similarity.shape)
            else:
                (i, j) = (0, 1)
            merged_id = "merged%d" % len(self)
            if i > j:
                (j,i) = (i,j)
            v = Motif.merge_observed(self[i].value, self[j].value, self.offset[i,j], self.antisense[i,j], merged_id)
            l = self[i]
            r = self[j]
            self.append(MotifTree(v, l, r))
            self.remove(l)
            self.remove(r)
            t2 = time.time()
            print "clustering... [%8.2f s] %d more trees left in MotifTreeList" % (t2 - t1, len(self))
        return self[0]
