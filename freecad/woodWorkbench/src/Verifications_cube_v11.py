# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay, delaunay_plot_2d
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.tri as mtri
from pylab import *
import os
import time
import numpy.matlib
import shutil
# from WoodMeshGenTools_v11 import *
import freecad.woodWorkbench.tools.WoodMeshGenTools_v11 as WoodMeshGen
from pathlib import Path
import FreeCAD as App

def verif():

    # performance check only
    startTime = time.time()

    # ==================================================================
    # Input parameters
    geoName = 'Debug_QML_v11_010824'
    outdir = App.ConfigGet('UserHomePath') + '/woodWorkbench'

    radial_growth_rule = 'binary'   
    iter_max = 10 # increase this number to achieve a more regular geometry
    print_interval = 500 # interval for printing prgress info

    # Radial cell growth parameters
    # length unit: mm
    r_min = 0   # inner radius of generation domain
    r_max = 50   # outer radius of generation domain
    nrings = 15 # number of rings
    width_heart = 0.15*(r_max-r_min)/nrings # ring width for the innermost ring
    width_early = 0.85*(r_max-r_min)/nrings # ring width for rings with early cells
    width_late = 0.15*(r_max-r_min)/nrings # ring width for rings with late cells
    generation_center = (0,0) # coordinates of generation domain center

    cellsize_early = 0.06
    cellsize_late = 0.03
    cellwallthickness_early = 0.01
    cellwallthickness_late = 0.006

    # clipping box parameters
    boundaryFlag = 'on'
    box_shape = 'square'
    # box_center = (1.005,0) # coordinates of clipping box center
    # box_size = 0.025 # side length

    box_center = (0.000,0) # coordinates of clipping box center
    box_size = 20 # side length

    # longitudinal direction parameters
    nsegments = 8 # increase this number to achieve a smoother transition in longitudinal direction if theta is large
    theta_min = 0 # unit: radian
    theta_max = 0.00 # unit: radian
    z_min = 0
    z_max = box_size # segment_length = (z_max - z_min) / nsegments
    long_connector_ratio = 0.2 # longitudinal joint length = ratio * segment_length

    # material parameters
    skeleton_density = 1.5e-9 # unit: tonne/mm3

    # generation parameters
    merge_operation = 'on'
    merge_tol = 0.015

    precrackFlag = 'off'
    precrack_widths = 0.1

    stlFlag = 'off'

    inpFlag = 'on'
    inpType = 'Abaqus'

    # ==================================================================
    # Remove directory/files if exists already
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

    if radial_growth_rule == 'binary':
        # ---------------------------------------------
        # binary radial growth pattern (e.g. wood microstructure with earlywood-latewood alternations)
        sites,radii = WoodMeshGen.CellPlacement_Binary(generation_center,r_max,r_min,nrings,width_heart,\
                            width_early,width_late,cellsize_early,cellsize_late,\
                            iter_max,print_interval)
    elif radial_growth_rule == 'binary_lloyd':
        # ---------------------------------------------
        # binary with Lloyd's algorithm (e.g. wood microstructure with earlywood-latewood alternations, but more regular cell shapes)
        sites,radii = WoodMeshGen.CellPlacement_Binary_Lloyd(geoName,outdir,generation_center,r_max,r_min,\
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
    # Clipping box (boundaries) of the model
    x_min,x_max,y_min,y_max,boundaries,boundary_points,boundarylines = \
        WoodMeshGen.Clipping_Box(box_shape,box_center,box_size,boundaryFlag)

    # if precracked
    x_indent_size = box_size*0.120
    y_indent_size = box_size*0.125
    x_precrack_size = box_size*0.1
    y_precrack_size = box_size*0.02
    x_indent = x_min + x_indent_size
    y_indent_min = box_center[1] - y_indent_size/2
    y_indent_max = box_center[1] + y_indent_size/2
    x_precrack = x_indent + x_precrack_size
    y_precrack = box_center[1]

    # Visualize the original 2D Voronoi diagram
    vor = Voronoi(sites[:,0:2])
    voronoi_plot_2d(vor, show_vertices=False,line_width=0.5, line_alpha=0.6, point_size=2)
    plt.xlim(x_min-0.1*abs(x_max-x_min), x_max+0.1*abs(x_max-x_min))
    plt.ylim(y_min-0.1*abs(y_max-y_min), y_max+0.1*abs(y_max-y_min))
    plt.plot(boundary_points[:, 0], boundary_points[:, 1], 'bo')
    # plt.show()

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

    # IGA beam parameters
    ninstance = 1
    nsvars_beam = 27 # number of svars for beam
    nsecgp = 4 # number of gauss points for beam sectional integration
    nsvars_secgp = 16 # number of svars at each gauss point
    iprops_beam = np.zeros(6)
    iprops_beam[0]    = 1 # section type index
    iprops_beam[1]    = 1 # number of instances
    iprops_beam[2]    = nconnector_t-nconnector_t_precrack # number of transverse connectors
    iprops_beam[3]    = nconnector_t_precrack # number of precracked transverse connectors
    iprops_beam[4]    = nconnector_l-nconnector_l_precrack # number of longitudinal connectors
    iprops_beam[5]    = nconnector_l_precrack # number of precracked longitudinal connectors
    iprops_beam = [int(x) for x in iprops_beam] 

    props_beam = np.zeros(16)
    props_beam[0]     = 1.5E-7 # Wood substance density [tonne/mm^3]
    props_beam[1]     = 0.8E+4 # Mesoscale elastic modulus [MPa]
    props_beam[2]     = 0.3 # Beam Poisson's ratio
    props_beam[3]     = cellwallthickness_early # Cross-sectional height [mm]
    props_beam[4]     = cellwallthickness_early # Cross-sectional width [mm]
    props_beam[5]     = 100 # Tensile Strength [MPa]
    props_beam[6]     = 200 # Tensile fracture energy [mJ/mm^2]
    props_beam[7]     = 4.1 # Shear Strength Ratio
    props_beam[8]     = 0.2 # Softening Exponent
    props_beam[9]     = 0.2 # Initial Friction
    props_beam[10]    = 0.0 # Asymptotic Friction
    props_beam[11]    = 600 # Transitional Stress [MPa]
    props_beam[12]    = 0.0 # Tensile Unloading
    props_beam[13]    = 0.0 # Shear Unloading
    props_beam[14]    = 0.0 # Shear Softening
    props_beam[15]    = 1.0 # Elastic Analysis Flag

    # Transverse connector parameters
    nsvars_conn_t = 32  # number of svars for transverse connector
    iprops_connector_t_bot = np.zeros(7)
    props_connector_t_bot = np.zeros(26)

    props_connector_t_bot[0]     = 1.5E-7 # Wood substance density [tonne/mm^3]
    props_connector_t_bot[1]     = 0.8E+4 # Mesoscale elastic modulus [MPa]
    props_connector_t_bot[2]     = 0.25 # Shear-Normal coupling coefficient
    props_connector_t_bot[3]     = height_connector_t # Connector height [mm]
    props_connector_t_bot[4]     = 0 # M-Distance [mm]
    props_connector_t_bot[5]     = -height_connector_t/2 # L-Distance [mm]
    props_connector_t_bot[6]     = 30.0 # Tensile Strength [MPa]
    props_connector_t_bot[7]     = 1.0 # Tensile characteristic length [mm] will be updated to # Tensile fracture energy [mJ/mm^2]
    props_connector_t_bot[8]     = 2.6 # Shear Strength Ratio
    props_connector_t_bot[9]     = 0.2 # Softening Exponent
    props_connector_t_bot[10]    = 0.2 # Initial Friction
    props_connector_t_bot[11]    = 0.0 # Asymptotic Friction
    props_connector_t_bot[12]    = 600 # Transitional Stress [MPa]
    props_connector_t_bot[12]    = 0.0 # Tensile Unloading
    props_connector_t_bot[14]    = 0.0 # Shear Unloading
    props_connector_t_bot[15]    = 0.0 # Shear Softening
    props_connector_t_bot[16]    = 0.0 # Elastic Analysis Flag
    props_connector_t_bot[17]    = 0.2 # Compressive Yielding Strength [MPa]
    props_connector_t_bot[18]    = 600 # Initial Hardening Modulus Ratio
    props_connector_t_bot[19]    = 0.0 # Transitional Strain Ratio
    props_connector_t_bot[20]    = 0.0 # Deviatoric Strain Threshold Ratio
    props_connector_t_bot[21]    = 0.0 # Deviatoric Damage Parameter
    props_connector_t_bot[22]    = 0.0 # Final Hardening Modulus Ratio
    props_connector_t_bot[23]    = 0.0 # Densification Ratio
    props_connector_t_bot[24]    = 0.0 # Volumetric Deviatoric Coupling
    props_connector_t_bot[25]    = 0.0 # Compressive Unloading

    iprops_connector_t_bot = [3,ninstance,nbeamElem,nconnector_t-nconnector_t_precrack,nconnector_t_precrack,nconnector_l-nconnector_l_precrack,nconnector_l_precrack]
    iprops_connector_t_bot = [int(x) for x in iprops_connector_t_bot] 
    iprops_connector_t_top = [2,ninstance,nbeamElem,nconnector_t-nconnector_t_precrack,nconnector_t_precrack,nconnector_l-nconnector_l_precrack,nconnector_l_precrack]
    iprops_connector_t_top = [int(x) for x in iprops_connector_t_top] 
    iprops_connector_t_reg = [1,ninstance,nbeamElem,nconnector_t-nconnector_t_precrack,nconnector_t_precrack,nconnector_l-nconnector_l_precrack,nconnector_l_precrack]
    iprops_connector_t_reg = [int(x) for x in iprops_connector_t_reg] 

    props_connector_t_reg = np.copy(props_connector_t_bot)
    props_connector_t_reg[3] = height_connector_t*2
    props_connector_t_reg[5] = 0
    props_connector_t_top = np.copy(props_connector_t_bot)
    props_connector_t_top[5] = height_connector_t/2

    # Longitudinal connector parameters
    nsvars_conn_l = 32  # number of svars for transverse connector
    iprops_connector_l = np.zeros(7)
    props_connector_l = np.zeros(24)

    props_connector_l[0]     = 1.5E-7 # Wood substance density [tonne/mm^3]
    props_connector_l[1]     = 0.8E+4 # Mesoscale elastic modulus [MPa]
    props_connector_l[2]     = 0.25 # Shear-Normal coupling coefficient
    props_connector_l[3]     = cellwallthickness_early*cellwallthickness_late # Connector sectional area [mm^2]
    props_connector_l[4]     = 20 # Tensile Strength [MPa]
    props_connector_l[5]     = long_connector_ratio*segment_length*1.05# 0.0105 # Tensile characteristic length [mm] will be updated to # Tensile fracture energy [mJ/mm^2]
    props_connector_l[6]     = 2.6 # Shear Strength Ratio
    props_connector_l[7]     = 0.2 # Softening Exponent
    props_connector_l[8]     = 0.2 # Initial Friction
    props_connector_l[9]     = 0.0 # Asymptotic Friction
    props_connector_l[10]    = 600 # Transitional Stress [MPa]
    props_connector_l[11]    = 0.0 # Tensile Unloading
    props_connector_l[12]    = 0.0 # Shear Unloading
    props_connector_l[13]    = 0.0 # Shear Softening
    props_connector_l[14]    = 1.0 # Elastic Analysis Flag
    props_connector_l[15]    = 0.2 # Compressive Yielding Strength [MPa]
    props_connector_l[16]    = 600 # Initial Hardening Modulus Ratio
    props_connector_l[17]    = 0.0 # Transitional Strain Ratio
    props_connector_l[18]    = 0.0 # Deviatoric Strain Threshold Ratio
    props_connector_l[19]    = 0.0 # Deviatoric Damage Parameter
    props_connector_l[20]    = 0.0 # Final Hardening Modulus Ratio
    props_connector_l[21]    = 0.0 # Densification Ratio
    props_connector_l[22]    = 0.0 # Volumetric Deviatoric Coupling
    props_connector_l[23]    = 0.0 # Compressive Unloading

    iprops_connector_l = [1,ninstance,nbeamElem,nconnector_t-nconnector_t_precrack,nconnector_t_precrack,nconnector_l-nconnector_l_precrack,nconnector_l_precrack]
    iprops_connector_l = [int(x) for x in iprops_connector_l] 

    # Strain rate effect parameters
    props_strainrate = np.zeros(4)
    props_strainrate[0] = 0.0 # Strain rate effect flag 
    props_strainrate[1] = 5.0 # Physical time scaling factor
    props_strainrate[2] = 1.0E-5 # Strain rate effect constant c0
    props_strainrate[3] = 5.0E-2 # Strain rate effect constant c1


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
    verif()
