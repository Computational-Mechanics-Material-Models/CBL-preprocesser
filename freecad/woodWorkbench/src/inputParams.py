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
## Primary Authors: Susan Alexis Brown, Hao Yin
## ===========================================================================

def inputParams(form):

    geoName = (form[0].geoName.text() or 'cube_debug')

    # Cell Growth Parameters
    radial_growth_rule = form[0].radial_growth_rule.currentText() #'binary'  #
    r_min = float(form[0].r_min.text() or 0)   # inner radius of generation domain
    r_max = float(form[0].r_max.text() or 1)   # outer radius of generation domain
    nrings = int(form[0].nrings.text() or 2) # number of rings
    width_heart = 0.15*(r_max-r_min)/nrings # ring width for the innermost ring
    width_early = 0.85*(r_max-r_min)/nrings # ring width for rings with early cells
    width_late = 0.15*(r_max-r_min)/nrings # ring width for rings with late cells
    generation_center = (0,0) # coordinates of generation domain center

    cellsize_early = float(form[0].cellsize_early.text() or 0.06)
    cellsize_late = float(form[0].cellsize_late.text() or 0.03)
    cellwallthickness_early = float(form[0].cellwallthickness_early.text() or 0.01)
    cellwallthickness_late = float(form[0].cellwallthickness_late.text() or 0.006)

    skeleton_density = float(form[0].skeleton_density.text() or 1.5e-9) # unit: tonne/mm3
    random_noise = float(form[0].random_noise.text() or 0.0) # unit: ?


    # Model Parameters     
    iter_max = int(form[0].iter_max.text() or 10) # increase this number to achieve a more regular geometry
    print_interval = 500 #int(form[0].print_interval.text() or 500) # interval for printing prgress info

    box_center = eval(form[0].box_center.text() or "(0.0,0.0)") # coordinates of clipping box center
    box_height = float(form[0].box_height.text() or 0.5) # specimen length
    box_width = float(form[0].box_width.text() or 0.5) # side length
    box_depth = float(form[0].box_depth.text() or 0.5) # side length
    box_size = box_width

    # longitudinal direction parameters
    nsegments = int(form[0].nsegments.text() or 2) 
    theta_min = float(form[0].theta_min.text() or 0) # unit: radian
    theta_max = 0.00 # unit: radian
    z_min = 0
    z_max = box_height #box_size # segment_length = (z_max - z_min) / nsegments
    long_connector_ratio = float(form[0].long_connector_ratio.text() or 0.05) # longitudinal joint length = ratio * segment_length

    NURBS_degree = 2

    # Knot parameters
    knotFlag = form[0].knotFlag.currentText()
    a1 = float(form[0].a1.text() or -0.1)
    a2 = float(form[0].a2.text() or 0.1)
    m1 = float(form[0].m1.text() or 0.05)
    m2 = float(form[0].m2.text() or 0.05)

    # Precrack Parameters
    # notch and precrack are always centered on leftmost edge
    x_notch_size = float(form[0].x_indent_size.text() or 0) # depth of notch box_size*0.120
    y_notch_size = float(form[0].y_indent_size.text() or 0) # width of notch box_size*0.125
    
    # Flags
    precrackFlag = form[0].precrackFlag.currentText()
    # precrack_widths = float(form[0].precrack_widths.text() or 0.1) # not used
    precrack_size = float(form[0].x_precrack_size.text() or 0) # depth of precrack box_size*0.1

    boundaryFlag = form[0].boundaryFlag.currentText()
    box_shape = form[0].box_shape.currentText()
    box_shape = box_shape.lower().replace(' ','_')

    merge_operation = form[0].merge_operation.currentText() 
    merge_tol = float(form[0].merge_tol.text() or 0.015)

    stlFlag = form[0].stlFlag.currentText()
    inpFlag = form[0].inpFlag.currentText()
    inpType = form[0].inpType.currentText()
    visFlag = form[0].visFlag.currentText()

    return geoName, radial_growth_rule, iter_max, print_interval, \
        r_min, r_max, nrings, width_heart, width_early, width_late, generation_center, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_size, box_height,\
        nsegments, theta_max, theta_min, z_max, z_min, long_connector_ratio, \
        x_notch_size, y_notch_size, precrack_size, \
        skeleton_density, merge_operation, merge_tol, precrackFlag, \
        stlFlag, inpFlag, inpType, random_noise, NURBS_degree, box_width, box_depth, visFlag, \
        knotFlag, m1, m2, a1, a2