#!/usr/bin/env python

import os
import sys
from distutils.core import setup, Extension

def main():
    setup(name="NPS",
          version="1.3.2",
          description="NPS - Nucleosome Positioning from Sequencing",
          author=['Yong Zhang', 'Hyunjin Shin'],
          author_email=['zy@jimmy.harvard.edu', 'shin@jimmy.harvard.edu'],
          packages=['NPS'],
          package_dir={'NPS' : 'lib'},
          scripts=['bin/SeqTagProcess.py', 'bin/SeqTag.py'],

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
          )

if __name__ == '__main__':
    main()
