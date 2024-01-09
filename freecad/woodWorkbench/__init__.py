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

import os
from .version import __version__

ICONPATH = os.path.join(os.path.dirname(__file__), "gui/icons")
GUIPATH = os.path.join(os.path.dirname(__file__), "gui/ui")
SRCPATH = os.path.join(os.path.dirname(__file__), "src")