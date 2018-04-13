#!/usr/bin/env python

# from distutils.core import setup, Extension
from setuptools import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs
import os
import codecs
import re

#Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(os.path.dirname(__file__), 'contour3d.py'),
                 encoding='utf8') as version_file:
    metadata = dict(re.findall(r"""__([a-z]+)__ = "([^"]+)""", version_file.read()))


setup(
      name='contour3d',
      version=metadata['version'],
      py_modules=['contour3d'],
      description='Contour 3D algorithm for grid data.',
      #long_description=README + '\n\n' +  CHANGES,
      classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        ],
      author='Masakazu Matsumoto',
      author_email='vitroid@gmail.com',
      url='https://github.com/vitroid/contour/',
      keywords=['contour','3D'],
      license='MIT',
      install_requires=['numpy', 'yaplotlib>=0.1'],
      entry_points = {
              'console_scripts': [
                  'contour3d = contour3d:main',
              ]
          }
)
