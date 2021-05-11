# -*- coding: utf-8 -*-
"""VTUinterface: a python PVD/VTU API"""

import os
import codecs
import re

from setuptools import setup, find_packages

# find __version__ ############################################################

def read(*parts):
    """Read file data."""
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, *parts), "r") as fp:
        return fp.read()

def find_version(*file_paths):
    """Find version without importing module."""
    version_file = read(*file_paths)
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

VERSION = find_version("VTUinterface", "_version.py")


###############################################################################

README = open("README.md").read()



setup(name="VTUinterface",
      version=VERSION,
      maintainer="Jörg Buchwald",
      maintainer_email="joerg_buchwald@ufz.de",
      author="Jörg Buchwald",
      author_email="joerg.buchwald@ufz.de",
      url="https://github.com/joergbuchwald/VTUinterface",
      long_description=README,
      long_description_content_type="text/markdown",
      classifiers=["Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Visualization",
          "Topic :: Scientific/Engineering :: Physics",
          "Topic :: Scientific/Engineering :: Mathematics",
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.8"],
      license="BSD-3 -  see LICENSE.txt",
      platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
      include_package_data=True,
      install_requires=["lxml", "vtk", "pandas", "scipy"],
      py_modules=["vtuIO"],
      package_dir={'': 'VTUinterface'})
