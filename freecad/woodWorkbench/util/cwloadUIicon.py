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
## ===========================================================================

# pyright: reportMissingImports=false
import os
from PySide import QtCore, QtGui
from freecad.woodWorkbench import ICONPATH


def cwloadUIicon(form,icon):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - form: The form that the icon will be assigned to
    - icon: The name of the icon file
    --------------------------------------------------------------------------
    ### Outputs ###
    - An icon that is assigned to the form
    --------------------------------------------------------------------------
    """

    # Assign the icon to the form
    form.setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, icon))))