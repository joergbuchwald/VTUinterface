# -*- coding: utf-8 -*-
"""VTUinterface: a python PVD/VTU API"""

from setuptools import setup, find_packages

setup(name="VTUinterface",
      version=0.63,
      maintainer="Jörg Buchwald",
      maintainer_email="joerg_buchwald@ufz.de",
      author="Jörg Buchwald",
      author_email="joerg.buchwald@ufz.de",
      url="https://github.com/joergbuchwald/VTUinterface",
      classifiers=["Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Visualization",
          "Topic :: Scientific/Engineering :: Physics",
          "Topic :: Scientific/Engineering :: Mathematics",
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.8"],
      license="MIT -  see LICENSE.txt",
      platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
      include_package_data=True,
      install_requires=["lxml", "vtk", "pandas", "scipy"],
      py_modules=["vtuIO"],
      package_dir={'': 'src'})
