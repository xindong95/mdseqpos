#!/usr/bin/env python

import os
import sys
import stat
from distutils.core import setup, Extension

#DEFAULT NUMPY_PATH e.g. /usr/local/lib/python2.6/site-packages/numpy
NUMPY_PATH = os.path.join(sys.prefix,"lib","python"+sys.version[:3],
                          'site-packages','numpy')

def check_pkg_dependencies():
    #CHECK for dependencies:
    try:
        import django
    except ImportError, e:
        sys.stderr.write("CRITICAL: DJANGO 1.1.1 or greater must be installed\n")
        sys.exit(1)

    try:
        import numpy
        NUMPY_PATH = numpy.__path__[0]
        #print NUMPY_PATH
    except ImportError, e:
        sys.stderr.write("CRITICAL: numpy 1.3 or greater must be installed\n")
        sys.exit(1)

def main():
    if not float(sys.version[:3])>=2.6:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.6! python 2.6.6 is recommended!\n")
        sys.exit(1)
    check_pkg_dependencies()
    setup(name="mdseqpos",
          version="0.578",
          description="Motif finding tools",
          author='Ying Lei, Len Taing, et al',
          author_email='lentaing@jimmy.harvard.edu',
          url='http://liulab.dfci.harvard.edu/ap',
          package_dir={'mdseqpos' : 'lib'},
          packages=['mdseqpos'],
          package_data={'mdseqpos': ['django/*.html', 'django/static/*',
                                     'weblogo/*.*', 'weblogo/LICENSE', 
                                     'weblogo/seqlogo', 'weblogo/README', 
                                     'weblogo/cache/*', 'weblogo/img/*', 
                                     'weblogo/release/*', 
                                     'img/*', 'tools/*',
                                     'database/*']},
          scripts=['bin/MDSeqPos.py', 'bin/MotifScan.py'],
          ext_modules = [
              Extension('mdseqpos._seq',
                        sources = ['lib/seqpos/seq.c',
                                   'lib/seqpos/genomescan_Markov0_3.c',
                                   ],
                        include_dirs = ['lib/seqpos',
                                        NUMPY_PATH + '/core/include'
                                        ]
                        ),
              Extension('mdseqpos._MDmod',
                        sources = ['lib/MDmodule/MDmod.c',
                                   'lib/MDmodule/keywd.c',
                                   ],
                        include_dirs = ['lib/seqpos',
                                        NUMPY_PATH + '/core/include'
                                        ]
                        ),
              ],
                    
          classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Artistic License',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Database',
        ],
          requires=['numpy (>=1.3.0)', 'Django (>=1.1.1)'],
      )

    #POST DISTUTILS touchups:
    #   1. fix the file permissions for weblogo/seqlogo
    #must fix the file permissions for weblogo/seqlogo
    try:
        import mdseqpos
        file_exe_perm = stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | \
                        stat.S_IROTH | stat.S_IXOTH
        seqlogo_path = os.path.join(mdseqpos.__path__[0], 'weblogo','seqlogo')
        twoBitToFa_path = os.path.join(mdseqpos.__path__[0], 'tools',
                                       'twoBitToFa')
        os.chmod(seqlogo_path, file_exe_perm)
        os.chmod(twoBitToFa_path, file_exe_perm)
    except ImportError, e:
        pass
    
if __name__ == '__main__':
    main()
