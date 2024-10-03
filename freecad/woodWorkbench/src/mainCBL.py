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

import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.tri as mtri
from pylab import *
import os
import time

import shutil
import tempfile
from pathlib import Path
import triangle as tr
from scipy.spatial import Voronoi, voronoi_plot_2d
import FreeCAD as App # type: ignore
import FreeCADGui as Gui # type: ignore
import feminout.importVTKResults as importVTK # type: ignore
import ObjectsFem # type: ignore
import FemGui # type: ignore
from PySide import QtCore, QtGui # type: ignore

from freecad.woodWorkbench.src.inputParams import inputParams
from freecad.woodWorkbench.src.outputParams import outputParams
from freecad.woodWorkbench.src.outputLog import outputLog
from freecad.woodWorkbench.src.readInput import readInput
from freecad.woodWorkbench.src.chronoInput import chronoInput
from freecad.woodWorkbench.src.chronoInput import chronoAux

import freecad.woodWorkbench.tools.WoodMeshGenTools_v11 as WoodMeshGen
from freecad.woodWorkbench.tools.rand_field_generator import RandomField


def main(self):

    # performance check only
    startTime = time.time()

    # ==================================================================
    self.form[1].progressBar.setValue(10) 
    self.form[1].statusWindow.setText("Status: Reading Parameters.") 
    # ==================================================================
    # Input parameters 
    [geoName, radial_growth_rule, iter_max, print_interval, \
        r_min, r_max, nrings, width_heart, width_early, width_late, generation_center, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_size, box_height, \
        nsegments, theta_max, theta_min, z_max, z_min, long_connector_ratio, \
        x_notch_size, y_notch_size, precrack_size, \
        skeleton_density, merge_operation, merge_tol, precrackFlag, \
        stlFlag, inpFlag, inpType, random_noise, NURBS_degree, box_width, box_depth, visFlag, \
        knotFlag, m1, m2, a1, a2, Uinf]\
            = inputParams(self.form)
    
    precrack_widths = 0 # for future use
    # random field parameters - SA to implement as input eventually
    random_field = 'off'
    RF_dimension = 1
    RF_dist_types= ["TruncatedGaussian"]
    RF_dist_params = [[1.,0.5,0]]
    RF_corr_l = 0.05
    RF_sampling_type="MC"

    # for future implementation - SA
    # grain_length = 5 # mm
    # height = z_max-z_min
    # nsegments = int(height/grain_length)
    # z_max = 100
    

    # ==================================================================
    self.form[1].progressBar.setValue(15) 
    self.form[1].statusWindow.setText("Status: Initializing Files.") 
    # ==================================================================

    # Make output directory if does not exist
    outDir = self.form[1].outputDir.text()

    i = 1
    geoNamenew = geoName
    while os.path.exists(Path(outDir + '/' + geoNamenew)):
        i +=1
        geoNamenew = geoName + str(i)
    geoName = geoNamenew
    # if not os.path.exists(Path(outDir + '/' + geoName)):
    os.makedirs(Path(outDir + '/' + geoName))


    # Prep naming variables
    meshName = geoName + "_mesh"
    analysisName = geoName + "_analysis"
    materialName = geoName + "_material"
    dataFilesName = geoName + '_dataFiles'
    visualFilesName = geoName + '_visualFiles'

    # Write input parameters to log file for repeated use
    outputLog(geoName, radial_growth_rule, iter_max, print_interval, \
        r_min, r_max, nrings, width_heart, width_early, width_late, generation_center, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_size, box_height, \
        nsegments, theta_max, theta_min, z_max, z_min, long_connector_ratio, \
        x_notch_size, y_notch_size, precrack_size, \
        skeleton_density, merge_operation, merge_tol, precrackFlag, \
        stlFlag, inpFlag, inpType, random_noise, NURBS_degree, box_width, box_depth, visFlag, \
        knotFlag, m1, m2, a1, a2, Uinf)
        

    # Make new freecad document and set view if does not exisit
    try:
        docGui.activeView().viewAxonometric()
    except:
        App.newDocument(geoName)
        docGui = Gui.activeDocument()
        docGui.activeView().viewAxonometric()
    Gui.runCommand('Std_PerspectiveCamera',1)

    # Generate new analysis object with associated material
    try:
        test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(analysisName)[0] != None)
    except:
        test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(analysisName) != [])
    if test == False:
        # Analysis
        analysis_object = ObjectsFem.makeAnalysis(App.ActiveDocument,analysisName)
        # Store Material
        material_object = ObjectsFem.makeMaterialSolid(App.ActiveDocument, materialName)
        mat = material_object.Material
        mat['Name'] = materialName
        material_object.Material = mat
        analysis_object.addObject(material_object)

    # ==================================================================
    self.form[1].progressBar.setValue(20) 
    self.form[1].statusWindow.setText("Status: Placing Cells.") 
    # ==================================================================
    # Place cells with a specific radial growth pattern

    # [sites,radii] = genSites(self.form)
    if radial_growth_rule == 'binary':
        # ---------------------------------------------
        # binary radial growth pattern (e.g. wood microstructure with earlywood-latewood alternations)
        sites,radii = WoodMeshGen.CellPlacement_Binary(generation_center,r_max,r_min,nrings,width_heart,\
                            width_early,width_late,cellsize_early,cellsize_late,\
                            iter_max,print_interval)
    elif radial_growth_rule == 'binary_lloyd':
        # ---------------------------------------------
        # binary with Lloyd's algorithm (e.g. wood microstructure with earlywood-latewood alternations, but more regular cell shapes)
        sites,radii = WoodMeshGen.CellPlacement_Binary_Lloyd(geoName,outDir,generation_center,r_max,r_min,\
                                                    nrings,width_heart,width_early,width_late,\
                                                    cellsize_early,cellsize_late,iter_max,\
                                                    print_interval)            
    elif radial_growth_rule == 'regular_hexagonal':
        # ----------------------------------
        # hexagonal honeycomb-like geometry
        sites,radii = WoodMeshGen.CellPlacement_Honeycomb(generation_center,r_max,r_min,nrings,\
                            box_center,box_size,width_heart,\
                            width_early,width_late,\
                            cellsize_early,cellsize_late,\
                            cellwallthickness_early,cellwallthickness_late,\
                            iter_max,print_interval)
    elif os.path.splitext(radial_growth_rule)[1] == '.npy':
        # ----------------------------------
        # load saved cell sites and radii data
        print('Loading saved sites')
        sites_path = Path(os.path.dirname(os.path.abspath(__file__)))
        sites,radii = WoodMeshGen.ReadSavedSites(sites_path,radial_growth_rule)
    else:
        print('Growth rule: {:s} is not supported for the current version, please check the README for more details.'.format(radial_growth_rule))
        # print('Now exiting...')
        # # exit()

    placementTime = time.time()
    nParticles = len(sites)
    print('{:d} particles/cells placed in {:.3f} seconds'.format(nParticles, (placementTime - startTime)))

    # ==================================================================
    self.form[1].progressBar.setValue(30) 
    self.form[1].statusWindow.setText("Status: Defining Boundaries.") 
    # ==================================================================
    # Clipping box (boundaries) of the model
    
    # x_min,x_max,y_min,y_max,boundaries,boundary_points_original,boundarylines = \
    #     WoodMeshGen.Clipping_Box(box_shape,box_center,box_size, box_width, box_depth,boundaryFlag,x_notch_size,y_notch_size)
    x_min,x_max,y_min,y_max,boundaries,boundary_points_original,boundarylines = \
        WoodMeshGen.Clipping_Box(box_shape,box_center,box_size,box_width,box_depth,boundaryFlag)

    sites_in = []
    for i in range(0,sites.shape[0]):
        if WoodMeshGen.check_isinside_boundbox2D(np.append(sites[i,:],0),x_min,x_max,y_min,y_max):
            sites_in.append(sites[i,:])

    # ---------------------------------------------
    # # Build flow mesh                    
    delaunay_vertices = np.concatenate((np.array(sites_in), boundary_points_original))
    # Conforming Delaunay
    conforming_delaunay = tr.triangulate({'vertices': delaunay_vertices}, 'pq0Dec')
    WoodMeshGen.BuildFlowMesh(outDir,geoName,conforming_delaunay,nsegments,long_connector_ratio,z_min,z_max)

    # num_interiorflowpts = np.shape(delaunay_vertices)[0]
    # flowpts_bound = np.concatenate(conforming_delaunay['vertices'][num_interiorflowpts:,:], boundary_points_original)

    # ---------------------------------------------
    # # Build mechanical mesh 
    vorsci = Voronoi(conforming_delaunay['vertices']) # different package to calculate voronoi region info

    ttvertices, ttedges, ttray_origins, ttray_directions = tr.voronoi(conforming_delaunay.get('vertices'))

    voronoiTime = time.time() 
    print('Original Voronoi tessellation generated in {:.3f} seconds'.format(voronoiTime - placementTime))


    # ==================================================================
    self.form[1].progressBar.setValue(40) 
    self.form[1].statusWindow.setText("Status: Clipping Mesh.") 
    # ==================================================================
    # Rebuild the Voronoi mesh
    # if merge_operation in ['on','On','Y','y','Yes','yes']:
    #     [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
    #     boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
    #     nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
    #         WoodMeshGen.RebuildVoronoi_ConformingDelaunay_merge(ttvertices,ttedges,ttray_origins,ttray_directions,\
    #                                                             boundaries,merge_tol,boundaryFlag)
        
    #     RebuildvorTime = time.time() 
    #     print('Voronoi tessellation rebuilt and merged in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))
    # else:
    #     [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
    #     boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
    #     nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
    #         WoodMeshGen.RebuildVoronoi_ConformingDelaunay(ttvertices,ttedges,ttray_origins,ttray_directions,\
    #                                                     boundaries,boundaryFlag)
    #     RebuildvorTime = time.time()

    # print(np.shape(ttvertices),np.shape(ttedges),np.shape(ttray_origins),np.shape(ttray_directions),np.shape(boundaries))

    [voronoi_vertices,boundary_points,finite_ridges_new,\
        boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
        nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
        WoodMeshGen.RebuildVoronoi_ConformingDelaunay_New(ttvertices,ttedges,ttray_origins,ttray_directions,\
                                                        boundaries,boundaryFlag)
    
    # print(np.shape(voronoi_vertices),np.shape(boundary_points),np.shape(finite_ridges_new),np.shape(boundary_ridges_new),np.shape(voronoi_ridges))
    # print((nvertex),(nvertices_in),(nfinite_ridge),(nboundary_ridge),(nboundary_pts),(nboundary_pts_featured),(nridge))

    

    RebuildvorTime = time.time()

    # ---------------------------------------------
    # # Visualize the mesh

    plt.figure()
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')

    # Flow
    vertsD = np.array(conforming_delaunay['vertices'])
    ax.triplot(vertsD[:, 0], vertsD[:, 1], conforming_delaunay['triangles'], 'bo-',markersize=3.,linewidth=0.5)

    # Main cells
    for beg, end in voronoi_ridges.astype(int):
        x0, y0 = voronoi_vertices[beg, :]
        x1, y1 = voronoi_vertices[end, :]
        ax.fill(
            [x0, x1],
            [y0, y1],
            facecolor='none',
            edgecolor='k',
            linewidth=1.0,
        )
        ax.plot(
            [x0, x1],
            [y0, y1],'ko',markersize=1.)
        
    # plt.show()
    
    plt.savefig(Path(outDir + '/' + geoName + '/' + geoName + '.png'))


    # ==================================================================
    self.form[1].progressBar.setValue(50) 
    self.form[1].statusWindow.setText("Status: Writing Vertices.") 
    # ================================================================== 
    [voronoi_vertices_3D,nvertices_3D,nlayers,segment_length,nctrlpt_per_elem,nctrlpt_per_beam,nconnector_t_per_beam,\
           nconnector_t_per_grain,theta,z_coord,npt_per_layer,npt_per_layer_normal,finite_ridges_3D,boundary_ridges_3D] = \
    WoodMeshGen.LayerOperation(NURBS_degree,nsegments,theta_min,theta_max,finite_ridges_new,boundary_ridges_new,nfinite_ridge,nboundary_ridge,\
                   z_min,z_max,long_connector_ratio,nvertices_in,nboundary_pts,\
                   voronoi_vertices,nvertex,voronoi_ridges,nridge,generation_center,random_noise,knotFlag, m1, m2, a1, a2, Uinf, box_center,box_depth)
    

    # Insert mid and quarter points on the Voronoi ridges (can be used as potential failure positions on cell walls)
    [all_pts_3D,npt_per_layer_vtk,all_pts_2D,all_ridges] = \
    WoodMeshGen.RidgeMidQuarterPts(voronoi_vertices_3D,nvertex,nvertices_in,voronoi_ridges,\
                    finite_ridges_new,boundary_ridges_new,finite_ridges_3D,boundary_ridges_3D,nfinite_ridge,\
                    nboundary_ridge,nboundary_pts,nlayers,voronoi_vertices)


    # ===============================================
    # Generate a file for vertices and ridges info
    [all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges] = \
        WoodMeshGen.VertexandRidgeinfo(all_pts_2D,all_ridges,\
                        npt_per_layer,npt_per_layer_normal,\
                        npt_per_layer_vtk,nridge,geoName,radii,generation_center,\
                        cellwallthickness_early,cellwallthickness_late,inpType)
    

    # ==================================================================
    self.form[1].progressBar.setValue(60) 
    self.form[1].statusWindow.setText("Status: Extruding Cells.") 
    # ==================================================================
    # Extrude in the parallel-to-grain (longitudinal) direction
    
    [IGAvertices,beam_connectivity_original,nbeam_total,\
    beam_connectivity,nbeamElem,\
    connector_t_bot_connectivity,connector_t_top_connectivity,\
    connector_t_reg_connectivity,connector_l_connectivity,\
    nconnector_t,nconnector_l,connector_l_vertex_dict] = \
    WoodMeshGen.GenerateBeamElement(voronoi_vertices_3D,nvertices_3D,NURBS_degree,nctrlpt_per_beam,nctrlpt_per_elem,nsegments,theta_min,theta_max,\
                        z_min,z_max,long_connector_ratio,npt_per_layer,voronoi_vertices,\
                        nvertex,voronoi_ridges,nridge,generation_center,all_vertices_2D,max_wings,\
                        flattened_all_vertices_2D,all_ridges,random_noise,nconnector_t_per_beam,nconnector_t_per_grain)

    BeamTime = time.time() 
    print('{:d} beam elements generated in {:.3f} seconds'.format(nbeamElem, (BeamTime - RebuildvorTime)))


    # ==================================================================
    # Calculate random field
    if random_field in ['on','On','Y','y','Yes','yes']:
        # RF = RandomField(RF_dimension,RF_dist_types,RF_dist_params,geoName,RF_corr_l,"text",[z_min,z_max],RF_sampling_type)
        # RF = RandomField(dimension=RF_dimension,"none",RF_dist_types,np.array([[1]]),RF_dist_params,geoName,RF_corr_l,"text",x_range=[x_min,x_max],RF_sampling_type,"binary")
        RF = RandomField(dimension = 1, readFromFolder = "none", dist_types=RF_dist_types, CC = np.array([[1]]), \
                         dist_params=RF_dist_params, name=geoName,corr_l =RF_corr_l, corr_f='exponential_square',\
                         x_range=[x_min,x_max], sampling_type = RF_sampling_type, filesavetype="binary")
        # RF = RandomField(dist_types=["TruncatedGaussian"], dist_params=[[1.,0.5,0]], name = "FIELD_1D_TruncNormal", corr_l = 0.05,filesavetype="text", x_range=[x_min,x_max],y_range=[y_min,y_max],z_range=[z_min,z_max], sampling_type="MC")
    else:
        RF = []

    # ===============================================
    # Insert precracks
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        
        x_notch = x_min + x_notch_size
        # y_notch_min = box_center[1] - y_notch_size/2
        # y_notch_max = box_center[1] + y_notch_size/2
        x_precrack = x_notch + precrack_size
        y_precrack = box_center[1]
        
        precrack_nodes = np.array([[x_notch, y_precrack, x_precrack, y_precrack]])
        [precrack_elem,nconnector_t_precrack,nconnector_l_precrack] = \
            WoodMeshGen.insert_precracks(all_pts_2D,all_ridges,nridge,npt_per_layer,\
                                    npt_per_layer_normal,npt_per_layer_vtk,\
                                    nlayers,precrack_nodes,precrack_size,\
                                    cellsize_early,nsegments)
    else:
        precrack_nodes = []
        precrack_elem = []
        nconnector_t_precrack = 0
        nconnector_l_precrack = 0

    # ==================================================================
    self.form[1].progressBar.setValue(70) 
    self.form[1].statusWindow.setText("Status: Calculating Mesh Info.") 
    # ==================================================================
    # Calculate mesh info
    [ConnMeshData,conn_l_tangents,height_connector_t] = \
        WoodMeshGen.ConnectorMeshFile(geoName,IGAvertices,connector_t_bot_connectivity,\
                    connector_t_reg_connectivity,connector_t_top_connectivity,\
                    connector_l_connectivity,all_vertices_2D,\
                    max_wings,flattened_all_vertices_2D,nsegments,segment_length,\
                    nctrlpt_per_beam,theta,nridge,connector_l_vertex_dict,\
                    random_field,RF)

    # ===============================================
    # Calculate model properties
    [mass,bulk_volume,bulk_density,porosity] = \
        WoodMeshGen.ModelInfo(box_shape,boundary_points,z_min,z_max,skeleton_density,ConnMeshData)

    # ===============================================
    # Bezier extraction 
    knotVec = WoodMeshGen.BezierExtraction(NURBS_degree,nbeam_total)
    npatch = beam_connectivity_original.shape[0]

    WoodMeshGen.BezierBeamFile(geoName,NURBS_degree,nctrlpt_per_beam,\
                    nconnector_t_per_beam,npatch,knotVec)

    # ==================================================================
    self.form[1].progressBar.setValue(80) 
    self.form[1].statusWindow.setText("Status: Generating Model Files.") 
    # ==================================================================
    # Generate input files for numerical simulations

    [nsvars_beam, nsecgp, nsvars_secgp, iprops_beam, props_beam, \
        nsvars_conn_t, iprops_connector_t_bot, props_connector_t_bot, iprops_connector_t_reg, props_connector_t_reg, \
        iprops_connector_t_top, props_connector_t_top, \
        nsvars_conn_l, iprops_connector_l, props_connector_l, props_strainrate] = outputParams(nconnector_l,nconnector_l_precrack,nconnector_t,nconnector_t_precrack,\
                 cellwallthickness_early,cellwallthickness_late,height_connector_t,nbeamElem,\
                    long_connector_ratio,segment_length)

    timestep     = 5.0E-9
    totaltime    = 5.0E-3
    z_min = np.min(z_coord)
    z_max = np.max(z_coord)
    # boundary_conditions = ['Hydrostatic']
    boundary_conditions = ['Bottom','Top','Left','Right','Front','Back']
    BC_velo_dof = 1 # 1-x, 2-y, 3-z
    BC_velo_value = box_size*0.03 # mm/s
    WoodMeshGen.AbaqusFile(geoName,NURBS_degree,npatch,nsegments,IGAvertices,beam_connectivity,\
                        connector_t_bot_connectivity,connector_t_reg_connectivity,\
                        connector_t_top_connectivity,connector_l_connectivity,\
                        segment_length,props_beam,iprops_beam,props_connector_t_bot,iprops_connector_t_bot,\
                        props_connector_t_reg,iprops_connector_t_reg,props_connector_t_top,\
                        iprops_connector_t_top,props_connector_l,iprops_connector_l,props_strainrate,\
                        timestep,totaltime,boundary_conditions,BC_velo_dof,BC_velo_value,\
                        x_max,x_min,y_max,y_min,z_max,z_min,boundaries,nsvars_beam,nsvars_conn_t,\
                        nsvars_conn_l,nsecgp,nsvars_secgp,cellwallthickness_early,merge_operation,merge_tol,\
                        precrackFlag,precrack_elem)
    chronoAux(geoName,IGAvertices,beam_connectivity,NURBS_degree,nctrlpt_per_beam,nconnector_t_per_beam,npatch,knotVec)
    # if inpFlag in ['on','On','Y','y','Yes','yes']:
    #     if inpType in ['abaqus','Abaqus','ABQ','abq','ABAQUS','Abq']:
    #         WoodMeshGen.AbaqusFile(geoName,NURBS_degree,npatch,nsegments,IGAvertices,beam_connectivity,\
    #                     connector_t_bot_connectivity,connector_t_reg_connectivity,\
    #                     connector_t_top_connectivity,connector_l_connectivity,\
    #                     segment_length,props_beam,iprops_beam,props_connector_t_bot,iprops_connector_t_bot,\
    #                     props_connector_t_reg,iprops_connector_t_reg,props_connector_t_top,\
    #                     iprops_connector_t_top,props_connector_l,iprops_connector_l,props_strainrate,\
    #                     timestep,totaltime,boundary_conditions,BC_velo_dof,BC_velo_value,\
    #                     x_max,x_min,y_max,y_min,z_max,z_min,boundaries,nsvars_beam,nsvars_conn_t,\
    #                     nsvars_conn_l,nsecgp,nsvars_secgp,cellwallthickness_early,merge_operation,merge_tol,\
    #                     precrackFlag,precrack_elem)
    #     elif inpType in ['Project Chrono','project chrono','chrono', 'Chrono']:
    #         # chronoInput(self.form)
    #         chronoAux(geoName,IGAvertices,beam_connectivity,NURBS_degree,nctrlpt_per_beam,nconnector_t_per_beam,npatch,knotVec)
    #     else:
    #         np.save(Path(outDir + '/' + geoName + '/' + geoName + '_sites.npy'),sites)
    #         np.save(Path(outDir + '/' + geoName + '/' + geoName + '_radii.npy'),radii)
    #         print('Input files type: {:s} is not supported for the current version, please check the README for more details.'.format(inpType))
    #         print('Generated cells and rings info has been saved.')
    # else:
    #     np.save(Path(outDir + '/' + geoName + '/' + geoName + '_sites.npy'),sites)
    #     np.save(Path(outDir + '/' + geoName + '/' + geoName + '_radii.npy'),radii)
    #     print('Generated cells and rings info has been saved.')
    
    FileTime = time.time() 
    print('Files generated in {:.3f} seconds'.format(FileTime - BeamTime))


    # ==============================================
    # Generate log file for the generation
    WoodMeshGen.LogFile(geoName,iter_max,r_min,r_max,nrings,width_heart,width_early,width_late,\
            generation_center,box_shape,box_center,box_size,x_min,x_max,y_min,y_max,
            cellsize_early,cellsize_late,cellwallthickness_early,cellwallthickness_late,\
            merge_operation,merge_tol,precrackFlag,precrack_widths,boundaryFlag,\
            nsegments,segment_length,theta_min,theta_max,z_min,z_max,long_connector_ratio,\
            NURBS_degree,nctrlpt_per_beam,nconnector_t_precrack,nconnector_l_precrack,\
            nParticles,nbeamElem,skeleton_density,mass,bulk_volume,bulk_density,porosity,\
            stlFlag,inpFlag,inpType,radial_growth_rule,\
            startTime,placementTime,voronoiTime,RebuildvorTime,BeamTime,FileTime)
    
    # ==================================================================
    self.form[1].progressBar.setValue(90) 
    self.form[1].statusWindow.setText("Status: Generating Vizualiation Files.") 
    # ==================================================================
    # Generate Paraview visulization files
    if visFlag in ['on','On','Y','y','Yes','yes']:
        
        # plt.show()
        
        WoodMeshGen.VisualizationFiles(geoName,NURBS_degree,nlayers,npt_per_layer_vtk,all_pts_3D,\
                       nsegments,nridge,voronoi_ridges,all_ridges,nvertex,\
                       nconnector_t,nconnector_l,nctrlpt_per_beam,ConnMeshData,\
                       conn_l_tangents,all_vertices_2D)  

        App.activeDocument().addObject('App::DocumentObjectGroup',visualFilesName)
        App.activeDocument().getObject(visualFilesName).Label = 'Visualization Files'

        importVTK.insert(str(outDir + '/' + geoName + '/' + geoName + '_beams.vtu'),App.ActiveDocument.Name)
        CBLbeamsVTU = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_beams')
        CBLbeamsVTU.Label = geoName + '_beams'
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(CBLbeamsVTU)
        CBLbeamsVTU.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + '/' + geoName + '/' + geoName + '_beams.vtu'))
        CBLbeamsVTU.Visibility = True

        importVTK.insert(str(outDir + '/' + geoName + '/' + geoName + '_conns.vtu'),App.ActiveDocument.Name)
        CBLconnsVTU = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_conns')
        CBLconnsVTU.Label = geoName + '_conns'
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(CBLconnsVTU)
        CBLconnsVTU.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + '/' + geoName + '/' + geoName + '_conns.vtu'))
        CBLconnsVTU.Visibility = False
        
        importVTK.insert(str(outDir + '/' + geoName + '/' + geoName + '_conns_vol.vtu'),App.ActiveDocument.Name)
        CBLconnsvVTU = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_conns_vol')
        CBLconnsvVTU.Label = geoName + '_conns_vol'
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(CBLconnsvVTU)
        CBLconnsvVTU.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + '/' + geoName + '/' + geoName + '_conns_vol.vtu'))
        CBLconnsvVTU.Visibility = False

        importVTK.insert(str(outDir + '/' + geoName + '/' + geoName + '_vertices.vtu'),App.ActiveDocument.Name)    
        CBLvertsVTU = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_vertices')
        CBLvertsVTU.Label = geoName + '_vertices'
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(CBLvertsVTU)
        CBLvertsVTU.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + '/' + geoName + '/' + geoName + '_vertices.vtu'))
        CBLvertsVTU.Visibility = False

        # =================================================
        # Generate 3D model files
        if stlFlag in ['on','On','Y','y','Yes','yes']:
            WoodMeshGen.StlModelFile(geoName)

        # ==================================================================    
        self.form[1].progressBar.setValue(100) 
        self.form[1].statusWindow.setText("Status: Complete.") 
        # ==================================================================

        # Switch to FEM GUI
        App.ActiveDocument.recompute()

        Gui.Control.closeDialog()
        # Gui.activateWorkbench("FemWorkbench")
        FemGui.setActiveAnalysis(App.activeDocument().getObject(analysisName))

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)
        Gui.runCommand('Std_PerspectiveCamera',1)
    
    else:

        # ==================================================================    
        self.form[1].progressBar.setValue(100) 
        self.form[1].statusWindow.setText("Status: Complete.") 
        # ==================================================================
        
        # Switch to FEM GUI
        App.ActiveDocument.recompute()
        Gui.Control.closeDialog()


if __name__ == '__main__':
    main()
