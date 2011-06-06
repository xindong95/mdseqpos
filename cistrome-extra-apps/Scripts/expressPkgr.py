#!/usr/bin/env python
"""
SYNOPSIS: given a set of arguments like this-
1:foo:GSMXXX 2:bar:GSMYYY 3:baz:GSMZZZ

This script will:
1. try to retrieve the .cel files associated with the GSMs
2. write a file called phenotype.txt (associating the Sample file w/ the key
   and the group:
   Sample    Key    Group
   a.cel     foo    1
   b.cel     bar    2
3. create a zip file packaging all of the .cel files and phenotype.txt
"""
import os
import sys
import time
import zipfile
import subprocess
import shutil
from ftplib import FTP,error_temp

try:
    import zlib
    compression = zipfile.ZIP_DEFLATED
except:
    compression = zipfile.ZIP_STORED

USAGE = """USAGE: expressPkgr 1:foo:GSMXXX 2:bar:GSMYYY 3:baz:GSMZZZ ...
WHERE: 
   1, 2, 3 - are the group ids
   foo, bar, baz - are the keys
   GSMXXX, GSMYYY, GSMZZZ - are the GSMIDs
"""
#------------------------------------------------------------------------------
#NOTE: THIS SECTION IS RIPPED FROM entrezutils/models.py

import urllib
from xml.dom.minidom import parseString
from xml.dom.minidom import Node
from xml.parsers.expat import ExpatError

#from django.db import models

import json

class XMLWrapper:
    """parses and stores an XML tree as follows:
    <elm0>
       <elm1 x=1 y=2>
          <elm2><elm3>foo</elm3></elm2>
       </elm1>
    <elm4>
       <elm5>bar</elm5>
    </elm4>
    <elm6/>
    </elm0>
    RETURNS: the following dictionary
    {_tagName="elm0", _attribs={},
    _children:[{_tagName="elm1", _attribs= {x:1, y:2},
                _children:[{_tagName="elm2", _attribs={},
                            _children:[{_tagName="elm3", _attribs={},
                                        _children:["foo"]}]
                            }]
                },
                {_tagName="elm4", _attribs={},
                 _children=[{_tagName="elm5", _attribs={},
                             _children:["bar"]}]
                },
                {_tagName="elm6", _attribs={}, _children:[]}
                ]}
    """
    def __init__(self, root_node):
        """given an XML DOM Tree, constructs the dictionary as described
        above"""
        self.root = self._parseDOM(root_node)

    def _parseDOM(self, node):
        """Constructs a dictionary that represents the XML node"""
        if node.nodeType == Node.TEXT_NODE:
            return node.nodeValue.strip()
        
        tmp = {'_tagName':node.tagName}
        if node.hasAttributes():
            attribs = {}
            for i in range(node.attributes.length):
                attribs[node.attributes.item(i).name] = \
                                                 node.attributes.item(i).value
            #print attribs
            tmp['_attribs'] = attribs
        if node.hasChildNodes():
            if len(node.childNodes) == 1 and \
               node.childNodes[0].nodeType == Node.TEXT_NODE and \
               node.childNodes[0].nodeValue.strip() != '':
                tmp['_value'] = node.childNodes[0].nodeValue.strip()
                tmp['_children'] = None
            else:
                children = filter(lambda x: x != '',
                                  [self._parseDOM(c) for c in node.childNodes])
                if len(children):
                    tmp['_children'] = children
                else:
                    tmp['_children'] = None
        else:
            tmp['_children'] = None

        return tmp
        
    def getElementsByTagName(self, tagName):
        """Returns an array of elements that match the tagName"""
        tmp = []
        def fn(node):
            if node['_tagName'] == tagName:
                tmp.append(node)
        treeWalker(self.root, fn)
        return tmp

def treeWalker(XMLWrapperNode, fn):
    """Walks XMLWrapper tree in BFS order, and applies the fn"""
    fn(XMLWrapperNode)
    if XMLWrapperNode['_children']:
        for x in XMLWrapperNode['_children']:
            treeWalker(x, fn)

def GeoQuery(accession):
    """We can access geo data over HTTP using the following url:
    http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=XXX&view=quick&form=xml&targ=self
    where XXX is the GEO accession string, e.g. GSE20852, GSM1137

    This call will return an xml file, which we will have to parse and store
    in the instance.

    So this is a data-holder class for the GEO xml files that are being made.
    OF course, all of this is hidden from the user, so all he does is:
    tmp = GEO('GSE20852')

    NOTE: this technically isn't a DJANGO model, but i've put it in this module
    b/c it's the best place for now
    """
    URL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&view=quick&form=xml&targ=self" % accession

    f = urllib.urlopen(URL)
    #print f.read()
    dom = parseString(f.read())
    f.close()
    tmp = XMLWrapper(dom.documentElement)
    #print tmp
    return tmp

# END ENTREZUTILS stuff
#------------------------------------------------------------------------------

def geoDownload(file_url):
    """Given a geo ftp file url, e.g. 
    ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM301nnn/GSM301351/GSM301351.CEL.gz
    
    tries to download the file and returns the fp to the file
    """
    parts = file_url.split("/")
    #Gives something like this:
    #[u'ftp:', u'', u'ftp.ncbi.nih.gov', u'pub', u'geo', u'DATA', u'supplementary', u'samples', u'GSM301nnn', u'GSM301351', u'GSM301351.CEL.gz']
    
    server = parts[2]
    path = parts[3:-1]
    filename = parts[-1]

    #GET the file!
    ftp = FTP(server)
    ftp.login()
    ftp.cwd("/".join(path))
    f = open(filename, "wb")
    ftp.retrbinary("RETR %s" % filename, f.write)
    f.close()
    ftp.close()
    
    return f

def handleArg(arg, pheno_file):
    """given an (string) arg in the form- GROUP:KEY:GSMID, e.g. 1:foo:GSM301351
    this fn will try to retrieve the .cel file from geo, and add the 
    entry to the pheno_file
    """
    (group, key, gsmid) = arg.split(":")
    
    #Get the Geo data
    geo = GeoQuery(gsmid)
    supplementary_data = geo.getElementsByTagName("Supplementary-Data")[0]
    file_url = supplementary_data['_value']
    
    f = geoDownload(file_url)
    
    #check if we need to gunzip the file
    if f.name.split(".")[-1] == "gz":
        ret_code = subprocess.call(["gunzip", f.name])
        pheno_file.write("%s\t%s\t%s\n" % \
                             (".".join(f.name.split(".")[:-1]), key, group))
    else:
        pheno_file.write("%s\t%s\t%s\n" % (f.name, key, group))
                         
def main():
    if len(sys.argv) < 2:
        print USAGE
        sys.exit(-1)

    args = sys.argv[1:]
    
    cwd = os.getcwd()
    working_dir = os.path.join(cwd, "tmp_%s" % time.time())
    os.mkdir(working_dir)
    os.chdir(working_dir)
    

    pheno_file = open("phenotype.txt", "w")
    pheno_file.write("Sample\tKey\tGroup\n")
    
    for a in args:
        handleArg(a, pheno_file)

    pheno_file.close()

    #now try to zip the contents    
    files = os.listdir(".") #get the files before creating the zip
    zip_file = zipfile.ZipFile("package.zip", "w")
    #print files
    for f in files:
        zip_file.write(f, compress_type=compression)
    zip_file.close()
    
    #move the zip_file
    shutil.move("package.zip", "..")

    #delete the temporary directory
    os.chdir(cwd)
    if os.path.exists(working_dir):
        shutil.rmtree(working_dir, ignore_errors=True)

if __name__ == "__main__":
    try:
        main()
    except error_temp as e:
        sys.stderr.write("FTP connection: %s\n",e.value)
        sys.exit()
    except ExpatError:
        sys.stderr.write("Can't find the GSM entry! GSM ID may be wrong.\n")
        sys.exit()
