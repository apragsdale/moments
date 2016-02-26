# Importing these adds a 'bdist_mpkg' option that allows building binary
# packages on OS X.
try:
    import setuptools
    import bdist_mpkg
except ImportError:
    pass

import os,sys

import numpy.distutils.core as core

#
# Microsoft Visual C++ only supports C up to the version iso9899:1990 (C89).
# gcc by default supports much more. To ensure MSVC++ compatibility when using
# gcc, we need to add extra compiler args. This code tries to ensure such
# arguments are added *only* when we're using gcc.
#
import numpy.distutils
compiler = numpy.distutils.ccompiler.get_default_compiler()
for arg in sys.argv:
    if arg.startswith('--compiler'):
        compiler = arg.split('=')[1]
if compiler in ['unix','mingw32','cygwin']:
    extra_compile_args = []
    # RNG: This seems to cause problems on some machines. To test for
    # compatibility with VC++, uncomment this line.
    #extra_compile_args = ['-std="iso9899:1990"', '-pedantic-errors']
else:
    extra_compile_args = []

numpy.distutils.core.setup(name='moments',
                           version='1.0.0',
                           author='Simon Gravel, Ryan Gutenkunst, Julien Jouganous',
                           author_email='simon.gravel@mcgill.ca',
                           url='http://simongravel.lab.mcgill.ca/Home.html',
                           packages=['moments'],
                           package_data = {'tests':['IM.fs']},
                           license='BSD'
                           )
