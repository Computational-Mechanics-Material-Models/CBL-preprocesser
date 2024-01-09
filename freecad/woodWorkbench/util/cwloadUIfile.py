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

import os
import FreeCADGui as Gui #type: ignore
from freecad.woodWorkbench import GUIPATH


def cwloadUIfile(filename):

        """
        Variables:
        --------------------------------------------------------------------------
        ### Inputs ###
        - filename: Name of the file to be loaded
        --------------------------------------------------------------------------
        ### Outputs ###
        - A module that is loaded into the left pane of the FreeCAD gui
        --------------------------------------------------------------------------
        """

        # Load the ui file
        ui_file_A = os.path.join(GUIPATH,filename)
        a = Gui.PySideUic.loadUi(ui_file_A)

        return a