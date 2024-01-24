# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pylab import *
import os
import time
import shutil
from pathlib import Path
import FreeCAD as App

from freecad.woodWorkbench.src.inputParams import inputParams
from freecad.woodWorkbench.src.genSites import genSites
from freecad.woodWorkbench.src.clipBox import clipBox
from freecad.woodWorkbench.src.outputParams import outputParams

import freecad.woodWorkbench.tools.WoodMeshGenTools_v11 as WoodMeshGen


def main():

    # performance check only
    startTime = time.time()
    outdir = App.ConfigGet('UserHomePath') + '/woodWorkbench'


    # ==================================================================
    # Input parameters
    [geoName, radial_growth_rule, iter_max, print_interval, \
        r_min, r_max, nrings, width_heart, width_early, width_late, generation_center, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_size, \
        nsegments, theta_max, theta_min, z_max, z_min, long_connector_ratio, \
        skeleton_density, merge_operation, merge_tol, precrack_widths, precrackFlag, \
        stlFlag, inpFlag, inpType] = inputParams()
    
    dumForm = [geoName, radial_growth_rule, iter_max, print_interval, \
        r_min, r_max, nrings, width_heart, width_early, width_late, generation_center, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_size, \
        nsegments, theta_max, theta_min, z_max, z_min, long_connector_ratio, \
        skeleton_density, merge_operation, merge_tol, precrack_widths, precrackFlag, \
        stlFlag, inpFlag, inpType]


    # ==================================================================
    # Remove directory/files if exists already
    try:
        shutil.rmtree(outdir + '/' + geoName)
    except:
        pass

    if not os.path.exists(Path(outdir + '/' + geoName)):
        os.makedirs(Path(outdir + '/' + geoName))
    
    # Set Path
    filpath = os.path.dirname(__file__)
    print(filpath)    


    # ==================================================================
    # Place cells with a specific radial growth pattern

    [sites,radii] = genSites(dumForm,outdir)

    placementTime = time.time()
    nParticles = len(sites)
    print('{:d} particles/cells placed in {:.3f} seconds'.format(nParticles, (placementTime - startTime)))


    # ==================================================================
    # Clipping box (boundaries) of the model
    
    [x_min, x_max, y_min, y_max, \
     x_indent, y_indent_min, y_indent_max, x_precrack, y_precrack, \
        vor, boundaries, boundarylines, boundary_points] = clipBox(box_shape,box_center,box_size,boundaryFlag,sites)

    voronoiTime = time.time() 
    print('Original Voronoi tessellation generated in {:.3f} seconds'.format(voronoiTime - placementTime))
    ax = plt.gca()
    boundarylines = patches.Polygon(boundary_points,closed=True,linewidth=2,edgecolor='k',facecolor='none')
    ax.set_aspect('equal', adjustable='box')
    ax.add_patch(boundarylines)


    # ==================================================================
    # Rebuild the Voronoi mesh
    if merge_operation in ['on','On','Y','y','Yes','yes']:
        [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
        boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
        nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
            WoodMeshGen.RebuildVoronoi_merge(vor,sites,boundaries,generation_center,x_min,x_max,y_min,y_max,box_center,box_shape,merge_tol,boundaryFlag)
        
        RebuildvorTime = time.time() 
        print('Voronoi tessellation rebuilt and merged in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))
    else:
        [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
        boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
        nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
            WoodMeshGen.RebuildVoronoi(vor,sites,boundaries,generation_center,x_min,x_max,y_min,y_max,box_center,box_shape,boundaryFlag)
        
        RebuildvorTime = time.time() 
        print('Voronoi tessellation rebuilt in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))

    # ===============================================
    # Insert mid and quarter points on the Voronoi ridges (can be used as potential failure positions on cell walls)
    [all_pts_2D,all_ridges,npt_per_layer,npt_per_layer_normal,npt_per_layer_vtk] = \
        WoodMeshGen.RidgeMidQuarterPts(voronoi_vertices,nvertex,nvertices_in,voronoi_ridges,\
                        finite_ridges_new,boundary_ridges_new,nfinite_ridge,\
                        nboundary_ridge,nboundary_pts,nboundary_pts_featured)

    # ==================================================================        
    # Generate a file for vertices and ridges info
    [all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges] = \
        WoodMeshGen.VertexandRidgeinfo(all_pts_2D,all_ridges,\
                        npt_per_layer,npt_per_layer_normal,\
                        npt_per_layer_vtk,nridge,geoName,radii,generation_center,\
                        cellwallthickness_early,cellwallthickness_late)
        
    ###############################################################################
    # Extrude in the parallel-to-grain (longitudinal) direction
    NURBS_degree = 2
    [IGAvertices,vertex_connectivity,beam_connectivity_original,nbeam_total,\
    beam_connectivity,nbeamElem,nctrlpt_per_beam,nlayers,segment_length,connector_t_connectivity,\
    connector_t_bot_connectivity,connector_t_top_connectivity,\
    connector_t_reg_connectivity,connector_l_connectivity,nconnector_t_per_beam,\
    nconnector_t_per_grain,nconnector_t,nconnector_l,nconnector_total,\
    theta,z_coord,connector_l_vertex_dict] = \
        WoodMeshGen.GenerateBeamElement(NURBS_degree,nsegments,\
                            theta_min,theta_max,z_min,z_max,\
                            long_connector_ratio,npt_per_layer,voronoi_vertices,\
                            nvertex,voronoi_ridges,nridge,generation_center,\
                            all_vertices_2D,max_wings,flattened_all_vertices_2D,all_ridges)

    BeamTime = time.time() 
    print('{:d} beam elements generated in {:.3f} seconds'.format(nbeamElem, (BeamTime - RebuildvorTime)))

    # ==================================================================
    # Insert precracks
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        precrack_nodes = np.array([[x_indent, y_precrack, x_precrack, y_precrack]])
        [precrack_elem,nconnector_t_precrack,nconnector_l_precrack] = \
            WoodMeshGen.insert_precracks(all_pts_2D,all_ridges,nridge,npt_per_layer,\
                                    npt_per_layer_normal,npt_per_layer_vtk,\
                                    nlayers,precrack_nodes,precrack_widths,\
                                    cellsize_early)
    else:
        precrack_nodes = []
        precrack_elem = []
        nconnector_t_precrack = 0
        nconnector_l_precrack = 0

    # ==================================================================
    # Calculate mesh info
    [ConnMeshData,conn_l_tangents,height_connector_t] = \
        WoodMeshGen.ConnectorMeshFile(geoName,IGAvertices,connector_t_bot_connectivity,\
                    connector_t_reg_connectivity,connector_t_top_connectivity,\
                    connector_l_connectivity,all_vertices_2D,\
                    max_wings,flattened_all_vertices_2D,nsegments,segment_length,\
                    nctrlpt_per_beam,theta,nridge,connector_l_vertex_dict)

    # ==================================================================
    # Calculate model properties
    [mass,bulk_volume,bulk_density,porosity] = \
        WoodMeshGen.ModelInfo(box_shape,boundary_points,z_min,z_max,skeleton_density,ConnMeshData)

    # ==================================================================
    # Bezier extraction 
    knotVec = WoodMeshGen.BezierExtraction(NURBS_degree,nbeam_total)
    npatch = beam_connectivity_original.shape[0]

    mkBezierBeamFile = WoodMeshGen.BezierBeamFile(geoName,NURBS_degree,nctrlpt_per_beam,\
                    nconnector_t_per_beam,npatch,knotVec)

    # ==================================================================
    # Generate Paraview visulization files
    WoodMeshGen.VisualizationFiles(geoName,NURBS_degree,nlayers,npt_per_layer_vtk,all_pts_2D,\
                    segment_length,theta,z_coord,nsegments,nridge,\
                    voronoi_ridges,generation_center,all_ridges,nvertex,nconnector_t,\
                    nconnector_l,nctrlpt_per_beam,ConnMeshData,conn_l_tangents,\
                    all_vertices_2D,max_wings,flattened_all_vertices_2D)
        
    plt.savefig(Path(outdir + '/' + geoName + '/' + geoName + '.png'))

    # ==================================================================
    # Generate 3D model files
    if stlFlag in ['on','On','Y','y','Yes','yes']:
        WoodMeshGen.StlModelFile(geoName)

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
    if inpFlag in ['on','On','Y','y','Yes','yes']:
        if inpType in ['abaqus','Abaqus','ABQ','abq','ABAQUS','Abq']:
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
        else:
            np.save(Path(outdir + geoName + '/' + geoName + '_sites.npy'),sites)
            np.save(Path(outdir + geoName + '/' + geoName + '_radii.npy'),radii)
            print('Input files type: {:s} is not supported for the current version, please check the README for more details.'.format(inpType))
            print('Generated cells and rings info has been saved.')
            print('Now exitting...')
            exit()
    
    FileTime = time.time() 
    print('Files generated in {:.3f} seconds'.format(FileTime - BeamTime))

    # ==================================================================
    # Generate log file for the generation
    WoodMeshGen.LogFile(geoName,iter_max,r_min,r_max,nrings,width_heart,width_early,width_late,\
            generation_center,box_shape,box_center,box_size,x_min,x_max,y_min,y_max,
            cellsize_early,cellsize_late,cellwallthickness_early,cellwallthickness_late,\
            merge_operation,merge_tol,precrackFlag,precrack_widths,boundaryFlag,\
            segment_length,theta_min,theta_max,z_min,z_max,long_connector_ratio,\
            NURBS_degree,nctrlpt_per_beam,nconnector_t_precrack,nconnector_l_precrack,\
            nParticles,nbeamElem,skeleton_density,mass,bulk_volume,bulk_density,porosity,\
            stlFlag,inpFlag,inpType,radial_growth_rule,\
            startTime,placementTime,voronoiTime,RebuildvorTime,BeamTime,FileTime)
    
if __name__ == '__main__':
    main()
