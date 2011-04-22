#!/usr/bin/env python
# this installation script was originally written by Giles Hall @ DFCI

import os

def check_build_tools():
    path = os.environ['PATH']
    for cmd in ('make', 'swig'):
        fq = None
        for dir in path.split( ':' ):
            testcmd = os.path.join( dir, cmd )
            try:
                os.stat( testcmd )
            except:
                continue
            fq = testcmd
        if not fq:
            print "Could not find necessary build command: %s" % cmd
            return False
    print "Found make and swig."
    import numpy
    numpy.ones(9,dtype=int) * 4
    return True

def build_wrappers():
    cwd = os.getcwd()
    os.chdir( 'src' )
    os.system( 'make all' )
    try:
        os.rename('AffyFileParser/extension/AffyFileParser.py', 'MAT/AffyFileParser.py')
        os.rename('gaussian/extension/gaussian.py', 'MAT/gaussian.py')
    except:
        pass
    os.chdir( cwd )

def clean_wrappers():
    cwd = os.getcwd()
    os.chdir( 'src' )
    os.system( 'make clean' )
    try:
        os.unlink('MAT/AffyFileParser.py')
        os.unlink('MAT/gaussian.py')
    except:
        pass
    os.chdir( cwd )

def call_distutils():
    from distutils.core import setup, Extension
    setup(
        name="MAT",
        version="3.07312009",
        description="A model-based algorithm to reliably detect regions enriched by transcription factor Chromatin ImmunoPrecipitation on Affymetrix tiling arrays.",

        long_description =
        '''
We propose a novel analysis algorithm, MAT, to reliably detect regions enriched by transcription factor Chromatin ImmunoPrecipitation (ChIP) on Affymetrix tiling arrays (chip). MAT models the baseline probe behavior by considering probe sequence and copy number on each array. The correlation between the baseline probe model estimates and the observed measurements can be as high as 0.72. MAT standardizes the probe value via the probe model, eliminating the need for sample normalization. A novel scoring function is applied to the standardized data to identify the ChIP-enriched regions, which allows robust p-value and false discovery rate calculations. MAT can detect ChIP-regions from a single ChIP sample, multiple ChIP samples, or multiple ChIP samples with controls with increasing accuracy. Based on the mock ChIP samples provided by the ENCODE consortium, MAT achieved 100% accuracy (0 false positive and 0 false negative) for the target detection of those spike-in plasmids, which are 2-, 4-, 8-, or 256-fold enriched compared with the genomic background. Quantitatively, MAT yielded a 0.95 correlation coefficient between the spike-in DNA concentration and the predicted score. Upon further analysis, MAT identified more than 70% of the true targets at 5% FDR cutoff from a single ChIP sample. This is a valuable feature for quickly testing the protocols and antibodies for ChIP-chip, and easily identifying ChIP-chip samples with questionable quality.
        ''',

        platforms = 'Linux/Unix/Cygwin',
        license = 'the Artistic License',
        author="Liu Lab",
        author_email='weili@jimmy.harvard.edu',
        url="http://chip.dfci.harvard.edu/~wli/MAT",
        updater_email='Madeleine_Lemieux@dfci.harvard.edu',
        options={'build_ext':{'swig_opts':'-c++'}},
        packages=["MAT"],
        package_dir = {
           'MAT':'src/MAT',
           'gaussian':'src/gaussian',
           'AffyFileParser':'src/AffyFileParser',
        },
        ext_package='MAT',
        scripts=['scripts/MAT','scripts/MatRep'],
        ext_modules = [
            Extension( '_gaussian',
                sources = ['src/gaussian/gsl/gauss.c',
                    'src/gaussian/gsl/gaussinv.c',
                    'src/gaussian/extension/gaussian_wrap.cpp'
                    ],
                include_dirs=['src/gaussian/gsl']
            ),
            Extension( '_AffyFileParser',
                sources = ['src/AffyFileParser/affymetrix/BARFileData.cpp',
                    'src/AffyFileParser/affymetrix/BARFileWriter.cpp',
                    'src/AffyFileParser/affymetrix/BEDFileData.cpp',
                    'src/AffyFileParser/affymetrix/BEDFileWriter.cpp',
                    'src/AffyFileParser/affymetrix/BPMAPFileData.cpp',
                    'src/AffyFileParser/affymetrix/BPMAPFileWriter.cpp',
                    'src/AffyFileParser/affymetrix/CDFFileData.cpp',
                    'src/AffyFileParser/affymetrix/CELFileData.cpp',
                    'src/AffyFileParser/affymetrix/CHPFileData.cpp',
                    'src/AffyFileParser/affymetrix/CHPFileWriter.cpp',
                    'src/AffyFileParser/affymetrix/FileIO.cpp',
                    'src/AffyFileParser/affymetrix/FileWriter.cpp',
                    'src/AffyFileParser/affymetrix/MSKFileData.cpp',
                    'src/AffyFileParser/affymetrix/MSKFileWriter.cpp',
                    'src/AffyFileParser/extension/AffyFileParser_wrap.cpp',
                    ],
                extra_compile_args = ['-w'],
                include_dirs=['src/AffyFileParser/affymetrix']
            ),
        ],
        requires=['numpy (>1.0.4)','rpy (>1.0.2)']
    )

if check_build_tools():
    clean_wrappers()
    build_wrappers()
    call_distutils()
    clean_wrappers()


