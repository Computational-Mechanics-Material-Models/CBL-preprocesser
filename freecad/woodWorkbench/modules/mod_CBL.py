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
import sys
import time
import tempfile
import numpy as np
from pathlib import Path
import functools
import math
import shutil

import FreeCADGui as Gui
import FreeCAD as App
import Part
import Part,PartGui
import Mesh
import MeshPartGui, FreeCADGui
import MeshPart
import Mesh, Part, PartGui
import MaterialEditor
import ObjectsFem
import FemGui
import Fem
import femmesh.femmesh2mesh
import Spreadsheet


from PySide import QtCore, QtGui

from freecad.woodWorkbench                                     import ICONPATH
from freecad.woodWorkbench                                     import GUIPATH
from freecad.woodWorkbench                                     import SRCPATH

# from freecad.woodWorkbench.src                                 import Verifications_cube_v11
from freecad.woodWorkbench.src                                 import mainCBL

from freecad.woodWorkbench.util.cwloadUIfile                    import cwloadUIfile
from freecad.woodWorkbench.util.cwloadUIicon                    import cwloadUIicon

# from freecad.woodWorkbench.output.mkChronoInput                  import mkChronoInput
# from freecad.woodWorkbench.output.mkAbaqusInput                  import mkAbaqusInput



class genWindow_CBL:
    def __init__(self):

        self.form = []

        # Load UI's for Side Panel
        self.form.append(cwloadUIfile("ui_CBL_cellParams.ui"))       
        self.form.append(cwloadUIfile("ui_CBL_generation.ui"))

        # Label, Load Icons, and Initialize Panels
        self.form[0].setWindowTitle("Model Properties")
        self.form[1].setWindowTitle("Model Generation")
        
        # cwloadUIicon(self.form[0],"ldpmOutput.svg")
        # cwloadUIicon(self.form[1],"PartDesign_AdditiveBox.svg")

        # # Set initial output directory
        self.form[1].outputDir.setText(str(Path(App.ConfigGet('UserHomePath') + '/woodWorkbench')))
        
        # Run generation for CBL
        QtCore.QObject.connect(self.form[1].generate, QtCore.SIGNAL("clicked()"), self.generation)

    def getStandardButtons(self):

        # Only show a close button
        # def accept() in no longer needed, since there is no OK button
        return int(QtGui.QDialogButtonBox.Close)


    def openDir(self):

        path = App.ConfigGet('UserHomePath')

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getExistingDirectory(None, "Open Directory",path,QtGui.QFileDialog.Option.ShowDirsOnly)         
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getExistingDirectory(None, "Open Directory",path,QtGui.QFileDialog.Option.ShowDirsOnly) 

        if OpenName == "":                                                            # if not selected then Abort process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            self.form[5].outputDir.setText(OpenName)

        return OpenName

    def generation(self):

        # Run file
        mainCBL.main(self)

        print('Fin.')
        # print("CModel package written to: " + outDir + outName)

        
# What to do when "Close" Button Clicked
    def reject(self):
        try:
            Gui.ActiveDocument.resetEdit()
            Gui.Control.closeDialog()
        except:
            Gui.Control.closeDialog()


class gen_CBL_Class():
    """My new command"""

    def GetResources(self):
        return {"Pixmap"  : os.path.join(ICONPATH, "ldpmOutput.svg"), # the name of a svg file available in the resources
                "MenuText": "CBL/CSL Simulation Inputs",
                "ToolTip" : "Generation of CBL simulation inputs"}

    def Activated(self):

        Gui.Control.showDialog(genWindow_CBL())

        return

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True

Gui.addCommand("mod_CBL", gen_CBL_Class())