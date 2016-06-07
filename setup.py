#!/usr/bin/env python

from distutils.core import setup

setup(name='Anchor',
      version='1.0',
      description='Anchor for STORM',
      author='Marcel Stefko',
      author_email='marcel.stefko@epfl.ch',
      url='https://github.com/MStefko/anchor-for-STORM',
      package_dir = {'': 'src'},
      py_modules = ['Anchor'],#'ImageStack', 'Anchor', 'beadData', 'STORM', 'tSTORMdata'],
      requires = ['scipy(>=0.17.0)','numpy(>=1.8.2)','libtiff(>=0.4.0)','pandas']
     )
