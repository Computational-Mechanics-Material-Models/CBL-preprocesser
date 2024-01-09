## ================================================================================
## CHRONO WOOD WORKBENCH
##
## Copyright (c) 2024 
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be found
## in the LICENSE file at the top level of the distribution
##
## ================================================================================
## Developed by Cusatis Computational Serivces, Inc.
## Primary Authors: Matthew Troemner, Susan Alexis Brown
## ================================================================================
##
## Description coming soon...
##
##
## ================================================================================

from setuptools import setup
import os

version_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 
                            "freecad", "woodWorkbench", "version.py")
with open(version_path) as fp:
    exec(fp.read())

setup(name='freecad.woodWorkbench',
      version=str(__version__),
      packages=['freecad',
                'freecad.woodWorkbench'],
      maintainer="sabrown",
      maintainer_email="susan.alexis.b@gmail.com",
      url="TBD",
      description="CBL Workbench",
      install_requires=['numpy','math','time'], 
      include_package_data=True)
