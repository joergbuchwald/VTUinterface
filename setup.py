# -*- coding: utf-8 -*-
"""VTUinterface: a python and julia PVD/VTU API"""

from setuptools import setup, find_packages

setup(name="VTUinterface",
      version=0.03,
      maintainer="Jörg Buchwald",
      maintainer_email="joerg_buchwald@ufz.de",
      author="Jörg Buchwald",
      author_email="joerg.buchwald@ufz.de",
      url="https://github.com/joergbuchwald/VTUinterface.jl",
      license="MIT -  see LICENSE.txt",
      platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
      include_package_data=True,
      install_requires=["lxml", "vtk"],
      py_modules=["vtuIO"])
