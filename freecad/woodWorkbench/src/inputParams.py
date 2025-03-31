## ================================================================================
## CHRONO WOOD WORKBENCH
##
## Copyright (c) 2025
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be found
## in the LICENSE file at the top level of the distribution
##
## ================================================================================
## Developed by Cusatis Computational Serivces, Inc.
## Primary Authors: Susan Alexis Brown, Hao Yin
## ===========================================================================
from math import sqrt, ceil
import os
from pathlib import Path
from freecad.woodWorkbench.src.outputLog import outputLog

def inputParams(form):


    # Cell Growth Parameters ------------------------
    radial_growth_rule = form[0].radial_growth_rule.currentText().lower().replace(' ','') #'binary'  
    species = form[0].species.currentText().lower().replace(' ','')  
    if species == 'norwayspruce':
        speciesShort = 'spruce'
        width_heart = 0.15 # ring width for the innermost ring
        ring_width = float(form[0].ring_width.text() or 2)
        late_ratio = float(form[0].ring_ratio.text() or 0.3)
        width_early = (1-late_ratio)*ring_width # ring width for rings with early cells 0.85 
        width_late = late_ratio*ring_width # ring width for rings with late cells 0.15
        cellsize_early = float(form[0].cellsize_early.text() or 0.032)
        cellsize_late = float(form[0].cellsize_late.text() or 0.017)
        cellwallthickness_early = float(form[0].cellwallthickness_early.text() or 0.0035)
        cellwallthickness_late = float(form[0].cellwallthickness_late.text() or 0.011)
        cell_length = float(form[0].cell_length.text() or 3)
    else: # default to generic
        speciesShort = 'generic'
        width_heart = 0.15 # ring width for the innermost ring
        ring_width = float(form[0].ring_width.text() or 2)
        late_ratio = float(form[0].ring_ratio.text() or 0.3)
        width_early = (1-late_ratio)*ring_width # ring width for rings with early cells 0.85 
        width_late = late_ratio*ring_width # ring width for rings with late cells 0.15
        cellsize_early = float(form[0].cellsize_early.text() or 0.03)
        cellsize_late = float(form[0].cellsize_late.text() or 0.02)
        cellwallthickness_early = float(form[0].cellwallthickness_early.text() or 0.003)
        cellwallthickness_late = float(form[0].cellwallthickness_late.text() or 0.006)
        cell_length = float(form[0].cell_length.text() or 10)
        
    
    # random field parameters
    randomFlag = form[0].randomFlag.currentText()
    randomParams = {}
    randomParams['RF_dist_types'] = ([form[0].dist_types.currentText()] or ["TruncatedGaussian"])
    randomParams['RF_dist_params'] = ([[float(i) for i in form[0].dist_params.text().split(',')]]) # or [[1.,0.05,0]]) # or doesn't work for empty after converting to float - convert later?
    randomParams['RF_corr_l'] = float(form[0].corr_l.text() or 0.1)
    randomParams['RF_sampling_type'] = (form[0].sampling_type.currentText() or "MC")


    # Geometry Parameters ------------------------
    box_shape = form[1].box_shape.currentText().lower().replace(' ','')
    box_center = eval(form[1].box_center.text() or "(5.75,0.0)") # coordinates of clipping box center
    generation_center = (0,0) # coordinates of generation domain center ??? - SA?
    x_notch_size = 0 # depth of notch
    y_notch_size = 0 # width of notch
    geoFile = ''
    if box_shape == 'cube':
        box_size = float(form[1].cube_size.text() or 0.5) # all side lengths
        box_height = box_size
        box_width = box_size
        box_depth = box_size
    elif box_shape == 'rectangle':
        box_height = float(form[1].box_height.text() or 0.5) # specimen length
        box_width = float(form[1].box_width.text() or 0.5) # x side length
        box_depth = float(form[1].box_depth.text() or 0.5) # y side length
        box_size = max(box_depth, box_width, box_depth)
    elif box_shape == 'notchedsquare':
        box_height = float(form[1].notch_height.text() or 0.5) # specimen length
        box_width = float(form[1].notch_width.text() or 0.5) # x side length
        box_depth = float(form[1].notch_depth.text() or 0.5) # y side length
        box_size = max(box_depth, box_width, box_depth)
        # notch and precrack are always centered on leftmost edge
        x_notch_size = float(form[1].x_indent_size.text() or box_width*0.1) # depth of notch
        y_notch_size = float(form[1].y_indent_size.text() or box_depth*0.1) # width of notch
    elif box_shape == 'input':
        geoFile = form[1].geoFile.text()
        box_height = float(form[1].geo_height.text() or 0.05)
        box_size = float(form[1].geo_size.text() or 2)
        box_width = box_size
        box_depth = box_size
    else:
        form[1].stackedWidget.setCurrentIndex(4)
        print('Geometry not recognized')
        return
    # patch for unkown bug that causes duplicate vertices for 1mm box 2/24/24   
    if box_width == box_depth and box_width == 1:
        box_width += 1e-10
        box_depth += 1e-10
    max_diag = sqrt(box_size**2*2) # diagonal length of square
    max_rad = sqrt(box_center[0]**2 + box_center[1]**2) # radius of box center
    r_max =  max_rad + max_diag + 0.1 # outer radius of generation domain with a buffer of 0.1
    nrings = ceil(r_max/ring_width) # number of rings to generate
    
    nsegments = max(ceil(box_height/cell_length*3),2) # number of segments in the longitudinal direction
    # at least 2 segments, with 3 per cell length
   
    precrackFlag = form[1].precrackFlag.currentText()
    precrack_start = eval(form[1].precrack_start.text() or "(-0.25,0)") # starting coords of precrack
    precrack_end = eval(form[1].precrack_end.text() or "(-0.2,0)") # end of precrack

    iter_max = int(form[1].iter_max.text() or 2) # lloyd iterations 
    NURBS_degree = 2
    # longitudinal direction parameters
    theta_min = float(form[1].theta_min.text() or 0) # unit: radian
    theta_max = 0.00 # unit: radian
    z_min = 0
    z_max = box_height #box_size 
    long_connector_ratio = float(form[1].long_connector_ratio.text() or 0.05) # longitudinal joint length = ratio * segment_length
    
    # knot parameters
    knotFlag = form[1].knotFlag.currentText()
    knotParams = {}
    knotParams['Uinf'] = float(form[1].knot_flow.text() or 1)
    knotParams['a1'] = float(form[1].a1.text() or 0.4)
    knotParams['a2'] = float(form[1].a2.text() or 0.1)
    knotParams['m1'] = float(form[1].m1.text() or 0.05)
    knotParams['m2'] = float(form[1].m2.text() or 0.05)


    # Model Parameters ------------------------     
    boundaryFlag = form[2].boundaryFlag.currentText()
    flowFlag = form[2].flowFlag.currentText()
    mergeFlag = form[2].merge_operation.currentText() 
    merge_tol = cellsize_early

    rayFlag = form[2].rayFlag.currentText()
    inpType = (form[2].inpType.currentText()).lower()
    visFlag = form[2].visFlag.currentText()

    if radial_growth_rule == 'debug': # override
        nsegments = int(2) 
        box_center = eval( "(0.0,0.0)") # coordinates of clipping box center
        box_height = 0.05 # specimen length
        box_width = 0.05 # side length
        box_depth = 0.05 # side length
        box_size = box_width
        z_min = 0
        z_max = box_height #box_size # segment_length = (z_max - z_min) / nsegments
        randomFlag = 'Off'
        x_notch_size = box_size*0.1 # depth of notch
        y_notch_size = box_size*0.1 # width of notch
        precrack_start = (-0.025,0)
        precrack_end = (-0.02,0)
    
    
    # Make output directory
    outDir = form[3].outputDir.text()
    geoName = (form[0].geoName.text() or (speciesShort + str(box_size).replace('.','_') + 'mm_' + box_shape))
    i = 1
    geoNamenew = geoName + str(i)
    while os.path.exists(Path(outDir + '/' + geoNamenew)):
        i +=1
        geoNamenew = geoName + str(i)
    geoName = geoNamenew
    # if not os.path.exists(Path(outDir + '/' + geoName)):
    os.makedirs(Path(outDir + '/' + geoName))

    # Write input parameters to log file for repeated use
    outputLog(geoName, radial_growth_rule, species, width_heart, ring_width, late_ratio, \
              cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
                cell_length, randomFlag, randomParams, box_shape, box_center, box_size, box_height, \
                  box_width, box_depth, x_notch_size, y_notch_size, precrackFlag, precrack_start, precrack_end, \
                    iter_max, theta_min, long_connector_ratio, knotFlag, knotParams, \
                      boundaryFlag,flowFlag, mergeFlag, rayFlag, inpType, visFlag, outDir,geoFile)

    return geoName, radial_growth_rule, iter_max, \
        nrings, width_heart, width_early, width_late, generation_center, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_size, box_height, \
        nsegments, theta_max, theta_min, z_max, z_min, long_connector_ratio, \
        x_notch_size, y_notch_size, precrack_start,precrack_end, \
        mergeFlag, merge_tol, precrackFlag, \
        inpType, randomFlag, randomParams, NURBS_degree, box_width, box_depth, visFlag, \
        knotFlag, knotParams, outDir,flowFlag, rayFlag
