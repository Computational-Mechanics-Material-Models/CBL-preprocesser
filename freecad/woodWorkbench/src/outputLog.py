import datetime
from pathlib import Path
import FreeCAD as App # type: ignore


def outputLog(geoName, radial_growth_rule, iter_max, \
        r_max, nrings, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_height, \
        nsegments, theta_min, long_connector_ratio, \
        x_notch_size, y_notch_size, precrack_size, \
        merge_operation, merge_tol, precrackFlag, \
        stlFlag, inpFlag, inpType, box_width, box_depth, visFlag, \
        knotFlag, knotParams, randomFlag, randomParams):

    # Generate log file
    logfile = open(Path(App.ConfigGet('UserHomePath') + '/woodWorkbench' + '/' + geoName + '/' + geoName + '-input.cwPar'),'w')        
    logfile.write("""
        // ================================================================================
        // CHRONO WORKBENCH - github.com/XXX
        //
        // Copyright (c) 2024 
        // All rights reserved. 
        //
        // Use of the code that generated this file is governed by a BSD-style license that
        // can be found in the LICENSE file at the top level of the distribution and at
        // github.com/XXX/LICENSE
        //
        // ================================================================================
        // CBL Workbench Parameter File
        // ================================================================================
        // 
        // CBL Workbench developed by Cusatis Computational Services
        //
        // ================================================================================
        \n\n""")
    
    logfile.write('radial_growth_rule= ' + str(radial_growth_rule) + '\n')
    logfile.write('r_max= ' + str(r_max) + '\n')
    logfile.write('nrings= ' + str(nrings) + '\n')

    logfile.write('cellsize_early= ' + str(cellsize_early) + '\n')
    logfile.write('cellsize_late= ' + str(cellsize_late) + '\n')
    logfile.write('cellwallthickness_early=' + str(cellwallthickness_early) + '\n')
    logfile.write('cellwallthickness_late= ' + str(cellwallthickness_late) + '\n')

    logfile.write('iter_max= ' + str(iter_max) + '\n')

    logfile.write('box_center= ' + str(box_center) + '\n')
    logfile.write('box_height= ' + str(box_height) + '\n')
    logfile.write('box_width=' + str(box_width) + '\n')
    logfile.write('box_depth= ' + str(box_depth) + '\n')

    logfile.write('nsegments= ' + str(nsegments) + '\n')
    logfile.write('theta_min= ' + str(theta_min) + '\n')
    logfile.write('long_connector_ratio= ' + str(long_connector_ratio) + '\n')

    logfile.write('knotFlag= ' + str(knotFlag) + '\n')
    logfile.write('Uinf= ' + str(knotParams.get('Uinf')) + '\n')
    logfile.write('a1= ' + str(knotParams.get('a1')) + '\n')
    logfile.write('a2= ' + str(knotParams.get('a2')) + '\n') 
    logfile.write('m1= ' + str(knotParams.get('m1')) + '\n')
    logfile.write('m2= ' + str(knotParams.get('m2')) + '\n') 

    logfile.write('randomFlag= ' + str(randomFlag) + '\n')
    logfile.write('dist_type= ' + str(randomParams.get('RF_dist_types')) + '\n')
    logfile.write('dist_params= ' + str(randomParams.get('RF_dist_params')) + '\n')
    logfile.write('corr_l= ' + str(randomParams.get('RF_corr_l')) + '\n') 
    logfile.write('sampling_type= ' + str(randomParams.get('RF_sampling_type')) + '\n')

    logfile.write('x_notch_size= ' + str(x_notch_size) + '\n')
    logfile.write('y_notch_size= ' + str(y_notch_size) + '\n')   

    logfile.write('precrackFlag= ' + str(precrackFlag) + '\n')
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        logfile.write('precrack_size= ' + str(precrack_size) + '\n')

    logfile.write('boundaryFlag= ' + str(boundaryFlag) + '\n')
    logfile.write('box_shape= ' + str(box_shape) + '\n')

    logfile.write('merge_operation= ' + str(merge_operation) + '\n')
    if merge_operation in ['on','On','Y','y','Yes','yes']: 
        logfile.write('merge_tol= ' + str(merge_tol) + '\n')
        

    logfile.write('stlFlag= ' + str(stlFlag) + '\n')
    logfile.write('inpFlag= ' + str(inpFlag) + '\n')
    logfile.write('inpType= ' + str(inpType) + '\n')
    logfile.write('visFlag= ' + str(visFlag) + '\n')

    logfile.write('\n')
    
    logfile.close()
