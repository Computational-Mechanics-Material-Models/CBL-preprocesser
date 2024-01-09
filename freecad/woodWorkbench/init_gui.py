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
import FreeCADGui as Gui #type: ignore
import FreeCAD as App #type: ignore
from PySide import QtGui #type: ignore
from FreeCADGui import Workbench #type: ignore

# Paths to Import
from freecad.woodWorkbench import ICONPATH
from freecad.woodWorkbench import GUIPATH
from freecad.woodWorkbench import SRCPATH

# # Chrono Scripts to Import
from freecad.woodWorkbench.modules import mod_CBL



class WoodWorkbench(Gui.Workbench):
    """
    class which gets initiated at startup of the gui
    """

    MenuText = "Wood Workbench"
    ToolTip = "A workbench for building CBL models for Project Chrono"
    Icon = os.path.join(ICONPATH, "PartDesign_Line.svg")
    toolbox = ["mod_CBL"] # a list of command names 


    def Initialize(self):
        """
        This function is called at the first activation of the workbench.
        here is the place to import all the commands
        """
        
        App.Console.PrintMessage("Switching to Wood Workbench\n")
        App.Console.PrintMessage("A workbench for building CBL models for Project Chrono and other solvers.\n")

        self.appendToolbar("Tools", self.toolbox) # creates a new toolbar with your commands
        self.appendMenu("Tools", self.toolbox) # creates a new menu
        self.appendMenu(["Wood Workbench"], self.toolbox) # appends a submenu to an existing menu



    def Activated(self):
        '''
        code which should be computed when a user switch to this workbench
        '''
        pass

    def Deactivated(self):
        '''
        code which should be computed when this workbench is deactivated
        '''
        pass

    def ContextMenu(self, recipient):
        """This function is executed whenever the user right-clicks on screen"""
        # "recipient" will be either "view" or "tree"
        self.appendContextMenu("My commands", self.list) # add commands to the context menu


    def GetClassName(self): 
        # This function is mandatory if this is a full Python workbench
        # This is not a template, the returned string should be exactly "Gui::PythonWorkbench"
        return "Gui::PythonWorkbench"



Gui.addWorkbench(WoodWorkbench())



