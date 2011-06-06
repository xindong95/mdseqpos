#!/usr/bin/env python

import os
import sys
from distutils.core import setup, Extension
#from Cython.Distutils import build_ext

def main():
    if not float(sys.version[:3])>=2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4! python 2.6.2 is recommended!\n")
        sys.exit(1)

    setup(name="cistrome-extra-apps",
          version="1.0",
          description="Cistrome extra applications",
          author='Tao (Foo) Liu, Jian Ma, Len Taing, Jacquie Wentz, and Wenbo Wang',
          author_email='taoliu@jimmy.harvard.edu',
          url='http://cistrome.dfci.harvard.edu/',
          package_dir={'CistromeAP' : '.'},
          packages=[
              'CistromeAP',
              'CistromeAP.taolib','CistromeAP.jianlib',
              'CistromeAP.taolib.CoreLib',
              'CistromeAP.taolib.CoreLib.DB','CistromeAP.taolib.CoreLib.FeatIO',
              'CistromeAP.taolib.CoreLib.BasicStat','CistromeAP.taolib.CoreLib.WWW',
              'CistromeAP.taolib.CoreLib.Parser','CistromeAP.taolib.CoreLib.SeqIO',
              'CistromeAP.taolib.CoreLib.BinKeeper','CistromeAP.taolib.CoreLib.Algorithm',
              ],

          scripts=[
              # Tao's scripts
              'Scripts/bed_correlation.py',
              'Scripts/conservation_plot.py',
              'Scripts/count_probes_in_peaks.py',
              'Scripts/count_probes_in_ranges.py',
              'Scripts/drawBED.py',
              'Scripts/fq2fa.py',
              'Scripts/naive_call_peaks.py',
              'Scripts/qc_chIP_peak.py',
              'Scripts/qc_chIP_whole.py',
              'Scripts/randPos',
              'Scripts/wig_call_peaks2',
              'Scripts/wig_call_peaks.py',
              'Scripts/wig_correlation_in_bed_file.py',
              'Scripts/wig_correlation.py',
              'Scripts/wig_extract_chrom.py',
              'Scripts/wiggle_reformat.py',
              'Scripts/xyz2image.py',
              # Jaquie's codes
              'Scripts/bedToWig2.py',
              'Scripts/wigToBed.py',
              'Scripts/wigLiftover.py',
              'Scripts/standardize_wig.py',
              'Scripts/venn_diagram.py',              
              # Jian's codes
              'Scripts/prof_sort.py',
              'Scripts/PCGA.py',
              'Scripts/heatmapr',
              # Len's codes
              'Scripts/expressPkgr.py',
              ],

          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: Database',
            ],
          requires=['PIL','Bio']
          )

if __name__ == '__main__':
    main()

