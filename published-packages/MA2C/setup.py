#!/usr/bin/env python
# Time-stamp: <2009-05-15 13:50:32 Tao Liu>

"""Description

Setup script for MA2C

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@status:  beta
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

import os
import sys
from distutils.core import setup, Extension

def main():
    if not float(sys.version[:3])>=2.4:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.4!\n")
        sys.exit(1)
    setup(name="pMA2C",
          version="1.1.3",
          description="Model Based Analysis for 2-color ChIP-chip",
          author='Xiaole Shirley Liu\'s Lab',
          author_email='ma2c.support@gmail.com',
          url='http://liulab.dfci.harvard.edu/MA2C/',
          package_dir={'MA2C' : 'lib'},
          packages=['MA2C'],
          scripts=['./bin/ma2c'],
          classifiers=[
            'Development Status :: 5 - productive',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            ],
          )
if __name__ == '__main__':
    main()
