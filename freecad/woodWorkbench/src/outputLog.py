from pathlib import Path
import FreeCAD as App # type: ignore
from datetime import datetime

def outputLog(geoName, radial_growth_rule, species, width_heart, ring_width, late_ratio, \
              cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
                cell_length, randomFlag, randomParams, box_shape, box_center, box_height, \
                  box_width, box_depth, x_notch_size, y_notch_size, precrackFlag, precrack_size, \
                    iter_max, theta_min, long_connector_ratio, knotFlag, knotParams, \
                      boundaryFlag,flowFlag, mergeFlag, rayFlag, inpType, visFlag, outDir):

    # Generate log file
    logfile = open(Path(outDir + '/' + geoName + '/' + geoName + '-input.cwPar'),'w')        
    logfile.write("""
        // ================================================================================
        // CHRONO CBL WORKBENCH - github.com/Cusatis-Computational-Services/CBL-preprocesser
        //
        // Copyright (c) 2025 
        // All rights reserved. 
        //
        // Use of the code that generated this file is governed by a BSD-style license that
        // can be found in the LICENSE file at the top level of the distribution and at
        // github.com/Cusatis-Computational-Services/CBL-preprocesser/LICENSE
        //
        // ================================================================================
        // CBL Workbench Parameter File
        // ================================================================================
        // 
        // CBL Workbench developed by Cusatis Computational Services, Inc.
        //
        // ================================================================================
        \n\n""")
    
  
    logfile.write('START' + '\n')
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    logfile.write('current_time= ' + current_time + '\n')

    logfile.write('geoName=' + str(geoName) + '\n')
    logfile.write('radial_growth_rule=' + str(radial_growth_rule) + '\n')

    # cell parameters
    logfile.write('species=' + str(species) + '\n')
    logfile.write('width_heart=' + str(width_heart) + '\n')
    logfile.write('ring_width=' + str(ring_width) + '\n')
    logfile.write('late_ratio=' + str(late_ratio) + '\n')
    logfile.write('cellsize_early=' + str(cellsize_early) + '\n')
    logfile.write('cellsize_late=' + str(cellsize_late) + '\n')
    logfile.write('cellwallthickness_early=' + str(cellwallthickness_early) + '\n')
    logfile.write('cellwallthickness_late=' + str(cellwallthickness_late) + '\n')
    logfile.write('cell_length=' + str(cell_length) + '\n')
                  
    # random field parameters
    logfile.write('randomFlag=' + str(randomFlag) + '\n')
    logfile.write('dist_type=' + str(randomParams.get('RF_dist_types')) + '\n')
    logfile.write('dist_params=' + str(randomParams.get('RF_dist_params')) + '\n')
    logfile.write('corr_l=' + str(randomParams.get('RF_corr_l')) + '\n') 
    logfile.write('sampling_type=' + str(randomParams.get('RF_sampling_type')) + '\n')
    
    # specimen geometry
    logfile.write('box_shape=' + str(box_shape) + '\n')
    logfile.write('box_center=' + str(box_center) + '\n')
    logfile.write('box_height=' + str(box_height) + '\n')
    logfile.write('box_width=' + str(box_width) + '\n')
    logfile.write('box_depth=' + str(box_depth) + '\n')
    logfile.write('x_notch_size=' + str(x_notch_size) + '\n')
    logfile.write('y_notch_size=' + str(y_notch_size) + '\n')

    logfile.write('precrackFlag=' + str(precrackFlag) + '\n')
    logfile.write('precrack_size=' + str(precrack_size) + '\n')

    logfile.write('iter_max=' + str(iter_max) + '\n')
    logfile.write('theta_min=' + str(theta_min) + '\n')
    logfile.write('long_connector_ratio=' + str(long_connector_ratio) + '\n')

    logfile.write('knotFlag=' + str(knotFlag) + '\n')
    logfile.write('Uinf=' + str(knotParams.get('Uinf')) + '\n')
    logfile.write('a1=' + str(knotParams.get('a1')) + '\n')
    logfile.write('a2=' + str(knotParams.get('a2')) + '\n') 
    logfile.write('m1=' + str(knotParams.get('m1')) + '\n')
    logfile.write('m2=' + str(knotParams.get('m2')) + '\n') 

    # model parameters
    logfile.write('boundaryFlag=' + str(boundaryFlag) + '\n')
    logfile.write('flowFlag=' + str(flowFlag) + '\n')
    logfile.write('mergeFlag=' + str(mergeFlag) + '\n')
    logfile.write('rayFlag=' + str(rayFlag) + '\n')
    logfile.write('inpType=' + str(inpType) + '\n')
    logfile.write('visFlag=' + str(visFlag) + '\n')
    
    logfile.write('outDir=' + str(outDir) + '\n')

    logfile.write('\n')
    
    logfile.close()
