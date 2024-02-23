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
from WoodMeshGenTools_v10 import *
# from hexalattice.hexalattice import *
from pathlib import Path


cwd = os.getcwd()
 
# performance check only
startTime = time.time()

# ==================================================================
# Input parameters

geoName = 'Cube_Notched_8x8mm_a60_QML_v10_012723'
geoName = 'test_notched'
path = 'meshes'

radial_growth_rule = 'wood_binary_lloyd'
iter_max = 1
print_interval = 9999999999

# length unit: mm
r_min = 0   # lower bound of radius
r_max = 17.5 # upper bound of radius
nrings = 15  # number of rings
width_heart = 0.15*(r_max-r_min)/nrings
width_early = 0.85*(r_max-r_min)/nrings
width_late = 0.15*(r_max-r_min)/nrings
log_center = (0,0) # coordinates of log center in global system of reference
# box_center = (5.70,-5.70) # coordinates of box center in global system of reference
box_center = (3.0*np.sqrt(3),-3.0) # coordinates of box center in global system of reference
# box_center = (3,0) # coordinates of box center in global system of reference
box_size = 2

box_center = (3,0) # coordinates of box center in global system of reference
box_size = 1

# if grooved
x_indent_size = box_size*10/38 # 10.0 #box_size*0.120
y_indent_size = box_size*6/38 #6.0 #box_size*0.125

# if precracked
x_precrack_size = 5/38*box_size
y_precrack_size = x_precrack_size*0.01

x_min = box_center[0] - box_size/2
x_max = box_center[0] + box_size/2
y_min = box_center[1] - box_size/2
y_max = box_center[1] + box_size/2
x_indent = x_min + x_indent_size
y_indent_min = box_center[1] - y_indent_size/2
y_indent_max = box_center[1] + y_indent_size/2
x_precrack = x_indent + x_precrack_size
y_precrack = box_center[1]


cellsize_early = 0.06
cellsize_late = 0.03
cellwallthickness_early = 0.01
cellwallthickness_late = 0.006

merge_operation = 'on'
merge_tol = 0.025

precrackFlag = 'on'
precrack_widths = 0.01

boundaryFlag = 'on'

stlFlag = 'on'
    
inpFlag = 'on'
inpType = 'Abaqus'
# ===============================================
# Boundaries of the model
box_shape = 'precracked_square'
# precracked
# ----------------------------------
boundary_points = np.array([[x_min, y_min],[x_max, y_min],[x_max, y_max],[x_min, y_max],[x_min, y_indent_max],[x_indent, y_indent_max],\
                            [x_indent,y_indent_min],[x_min,y_indent_min]])

l0 = [(x_min, y_min), (x_max, y_min)]
l1 = [(x_max, y_min), (x_max, y_max)]
l2 = [(x_max, y_max), (x_min, y_max)]
l3 = [(x_min, y_max), (x_min, y_indent_max)]
l4 = [(x_min, y_indent_max), (x_indent, y_indent_max)]
l5 = [(x_indent, y_indent_max), (x_indent,y_indent_min)]
l6 = [(x_indent,y_indent_min), (x_min,y_indent_min)]
l7 = [(x_min,y_indent_min), (x_min, y_min)]

boundaries = [('bottom',l0),('right',l1),('top',l2),('l3',l3),('l4',l4),('l5',l5),('l6',l6),('l7',l7)]

if boundaryFlag in ['on','On','Y','y','Yes','yes']:
    boundarylines = patches.Rectangle(boundaries[0][1][0], x_max-x_min, y_max-y_min, linewidth=2, edgecolor='k',facecolor='none')
else:
    boundarylines = patches.Rectangle(boundaries[0][1][0], x_max-x_min, y_max-y_min, linewidth=2, linestyle='-.', edgecolor='k',facecolor='none')

# ===============================================
# Remove directory/files if exists already
try:
    shutil.rmtree(Path(path, geoName))
except:
    pass

if not os.path.exists(Path(path, geoName)):
    os.makedirs(Path(path, geoName))
    
# ===============================================
# Generate annulus structure for wood microstructure
omega = 1.0
# sites,radius = Annualringseeds_lloyd(geoName,path,log_center,r_max,r_min,nrings,width_heart,width_early,width_late,\
#                   cellsize_early,cellsize_late,iter_max,print_interval,omega)

# user may save the generated "sites" and "radius" with a seed for repeated use (TBD)    
# sites = np.load('test_grooved_sites.npy')
# radius = np.load('test_grooved_radius.npy')
# sites = np.load('1mm_cube_I50000_v8_sites.npy')
# radius = np.load('1mm_cube_I50000_v8_radius.npy')
# sites = np.load('test_sites.npy')
# radius = np.load('test_radius.npy')
# sites = np.load('1x1_notched_I50000_v71_sites.npy')
# radius = np.load('1x1_notched_I50000_v71_radius.npy')
# sites = np.load('cube_1x1mm_v7_sites.npy')
# radius = np.load('cube_1x1mm_v7_radius.npy')
# sites = np.load('1mm_cube_I50000_220816_sites.npy')
# radius = np.load('1mm_cube_I50000_220816_radius.npy')
# sites = np.load('r8_n8_h02_e07_l02_ce0016_cl0008_I50000_sites.npy')
# radius = np.load('r8_n8_h02_e07_l02_ce0016_cl0008_I50000_radius.npy')
sites = np.load('r17p5_n15_h015_e085_l015_ce006_cl003_Iloyd_sites.npy')
radius = np.load('r17p5_n15_h015_e085_l015_ce006_cl003_Iloyd_radius.npy')

placementTime = time.time()
nParticles = len(sites)

print('{:d} particles/cells placed in {:.3f} seconds'.format(nParticles, (placementTime - startTime)))
                
# Visualize the original Voronoi diagram generated with the wood mesh sites    
vor = Voronoi(sites[:,0:2])
voronoi_plot_2d(vor, show_vertices=False,line_width=0.5, line_alpha=0.6, point_size=2)
# #------------------------------------------------------------------------------
# # Plot wood cell circles
# plt.plot(sites[:, 0], sites[:, 1], 'bo')
# circle_out = circles(sites[:, 0], sites[:, 1], sites[:, 2],alpha=0.5,edgecolor='none')
# colorbar(circle_out)
# #------------------------------------------------------------------------------
plt.xlim(x_min-0.1*abs(x_max-x_min), x_max+0.1*abs(x_max-x_min))
plt.ylim(y_min-0.1*abs(y_max-y_min), y_max+0.1*abs(y_max-y_min))
# plt.plot(boundary_points[:, 0], boundary_points[:, 1], 'bo')
plt.show()


voronoiTime = time.time() 
print('Original Voronoi tessellation generated in {:.3f} seconds'.format(voronoiTime - placementTime))

ax = plt.gca()
boundarylines = patches.Polygon(boundary_points,closed=True,linewidth=2,edgecolor='k',facecolor='none')
ax.set_aspect('equal', adjustable='box')
ax.add_patch(boundarylines)
# ===============================================
# Rebuild the Voronoi mesh
if merge_operation in ['on','On','Y','y','Yes','yes']:
    [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
     boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
     nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
        RebuildVoronoi_new(vor,sites,boundaries,log_center,x_min,x_max,y_min,y_max,box_center,merge_tol,boundaryFlag)
    
    RebuildvorTime = time.time() 
    print('Voronoi tessellation rebuilt and merged in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))
else:
    [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
      boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
      nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
        RebuildVoronoi(vor,sites,boundaries,log_center,x_min,x_max,y_min,y_max,box_center,boundaryFlag)
    
    RebuildvorTime = time.time() 
    print('Voronoi tessellation rebuilt in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))


# ===============================================
# The dual Delaunay the Voronoi mesh
# Visualize the original Delaunay triangulation generated with the wood mesh sites    
# tri = Delaunay(voronoi_vertices[:,0:2])
# triang = mtri.Triangulation(sites[:,0], sites[:,1])
# plt.triplot(triang, 'go-')
# delaunay_plot_2d(tri)


# plt.triplot(np.copy(sites[:,0]).flatten() ,np.copy(sites[:,1]).flatten(), tri.simplices)

# ===============================================
# Insert Mid and Quarter points on the Voronoi ridges (potential failure positions on cell walls)
[all_pts_2D,all_ridges,npt_per_layer,npt_per_layer_normal,npt_per_layer_vtk] = \
    RidgeMidQuarterPts(voronoi_vertices,nvertex,nvertices_in,voronoi_ridges,\
                       finite_ridges_new,boundary_ridges_new,nfinite_ridge,\
                       nboundary_ridge,nboundary_pts,nboundary_pts_featured)
        
# Generate a file for the geometry info for vertices and ridges
[all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges] = \
    VertexandRidgeinfo(all_pts_2D,all_ridges,\
                       npt_per_layer,npt_per_layer_normal,\
                       npt_per_layer_vtk,nridge,geoName,radius,log_center,\
                       cellwallthickness_early,cellwallthickness_late)
    
###############################################################################
# Extrusion in the parallel-to-grain (longitudinal) direction
fiberlength = 0.5*2.5#0.5*box_size
theta_min = 0 
theta_max = 0 #0.05*2*math.pi
z_min = 0
z_max = 2.5#box_size

long_connector_ratio = 0.02

NURBS_degree = 2
nctrlpt_per_beam = 5

[woodIGAvertices,vertex_connectivity,beam_connectivity_original,nbeam_total,\
 beam_connectivity,nbeamElem,nlayer,connector_t_connectivity,\
 connector_t_bot_connectivity,connector_t_top_connectivity,\
 connector_t_reg_connectivity,connector_l_connectivity,nconnector_t_per_beam,\
 nconnector_t_per_grain,nconnector_t,nconnector_l,nconnector_total,\
 theta,z_coord,nbeam_per_grain,connector_l_vertex_dict] = \
    GenerateBeamElement(NURBS_degree,nctrlpt_per_beam,\
                        fiberlength,theta_min,theta_max,z_min,z_max,\
                        long_connector_ratio,npt_per_layer,voronoi_vertices,\
                        nvertex,voronoi_ridges,nridge,log_center,\
                        all_vertices_2D,max_wings,flattened_all_vertices_2D,all_ridges)

BeamTime = time.time() 
print('{:d} beam elements generated in {:.3f} seconds'.format(nbeamElem, (BeamTime - RebuildvorTime)))

########################### Insertion of precracks ############################
if precrackFlag in ['on','On','Y','y','Yes','yes']:
    # pre-crack node list: [c1x1 c1y1 c1x2 c1y2 ...]
    precrack_nodes = np.array([[x_indent, y_precrack, x_precrack, y_precrack]])
    precrack_elem = insert_precracks(all_pts_2D,all_ridges,nridge,npt_per_layer,\
                                       npt_per_layer_normal,npt_per_layer_vtk,\
                                       nlayer,precrack_nodes,precrack_widths,\
                                       cellsize_early,nbeam_per_grain)
else:
    precrack_nodes = []
    precrack_elem = []

nconnector_t_precrack = len(precrack_elem)*3 # times 3 for bot,reg,top layers
nconnector_l_precrack = 0

######################### Connector Calculations ##############################
height_connector_t = fiberlength/4 # fiberlength/((nconnector_t_per_grain-1)*2)

ConnMeshData = ConnectorMeshFile(geoName,woodIGAvertices,connector_t_bot_connectivity,\
                 connector_t_reg_connectivity,connector_t_top_connectivity,\
                 height_connector_t,connector_l_connectivity,all_vertices_2D,\
                 max_wings,flattened_all_vertices_2D,nbeam_per_grain,nridge,\
                 connector_l_vertex_dict)
    
########################### Bezier extraction #################################
knotVec = BezierExtraction(NURBS_degree,nctrlpt_per_beam,nbeam_total)
npatch = beam_connectivity_original.shape[0]

mkBezierBeamFile = BezierBeamFile(geoName,NURBS_degree,nctrlpt_per_beam,\
                   nconnector_t_per_beam,npatch,knotVec)

########################### Paraview VTK data #################################
VisualizationFiles(geoName,NURBS_degree,nlayer,npt_per_layer_vtk,all_pts_2D,\
                   fiberlength,theta,z_coord,nbeam_per_grain,nridge,\
                   voronoi_ridges,log_center,all_ridges,nvertex,nconnector_t,\
                   nconnector_l,nctrlpt_per_beam,ConnMeshData,all_vertices_2D,\
                   max_wings,flattened_all_vertices_2D,precrack_nodes,precrack_elem,\
                   nconnector_t_precrack,nconnector_l_precrack)
    
############################# Abaqus inp File #################################

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
props_beam[0]     = 1.5E-9 # Wood substance density [tonne/mm^3]
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

props_connector_t_bot[0]     = 1.5E-9 # Wood substance density [tonne/mm^3]
props_connector_t_bot[1]     = 0.8E+4 # Mesoscale elastic modulus [MPa]
props_connector_t_bot[2]     = 0.25 # Shear-Normal coupling coefficient
props_connector_t_bot[3]     = height_connector_t # Connector height [mm]
props_connector_t_bot[4]     = 0 # M-Distance [mm]
props_connector_t_bot[5]     = -height_connector_t/2 # L-Distance [mm]
props_connector_t_bot[6]     = 30.0 # Tensile Strength [MPa]
props_connector_t_bot[7]     = 100.0 # Tensile characteristic length [mm] will be updated to # Tensile fracture energy [mJ/mm^2]
props_connector_t_bot[8]     = 2.6 # Shear Strength Ratio
props_connector_t_bot[9]     = 0.2 # Softening Exponent
props_connector_t_bot[10]    = 0.2 # Initial Friction
props_connector_t_bot[11]    = 0.2 # Asymptotic Friction
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
# height_connector_l = fiberlength/((nconnector_l_per_grain-1)*2)

props_connector_l[0]     = 1.5E-9 # Wood substance density [tonne/mm^3]
props_connector_l[1]     = 0.8E+4 # Mesoscale elastic modulus [MPa]
props_connector_l[2]     = 0.25 # Shear-Normal coupling coefficient
props_connector_l[3]     = cellwallthickness_early*cellwallthickness_late # Connector sectional area [mm^2]
props_connector_l[4]     = 20 # Tensile Strength [MPa]
props_connector_l[5]     = long_connector_ratio*fiberlength*1.05# 0.0105 # Tensile characteristic length [mm] will be updated to # Tensile fracture energy [mJ/mm^2]
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


timestep     = 2.0E-8
totaltime    = 3.0E-3
z_min = np.min(z_coord)
z_max = np.max(z_coord)
# boundary_conditions = ['Hydrostatic']
boundary_conditions = ['Bottom','Top','Left','Right','Front','Back']
BC_velo_dof = 2 # 1-x, 2-y, 3-z
BC_velo_value = box_size*1.25 # mm/s
    
mkAbaqusFile = AbaqusFile(geoName,NURBS_degree,npatch,nbeam_per_grain,woodIGAvertices,beam_connectivity,\
               connector_t_bot_connectivity,connector_t_reg_connectivity,\
               connector_t_top_connectivity,connector_l_connectivity,\
               fiberlength,props_beam,iprops_beam,props_connector_t_bot,iprops_connector_t_bot,\
               props_connector_t_reg,iprops_connector_t_reg,props_connector_t_top,\
               iprops_connector_t_top,props_connector_l,iprops_connector_l,props_strainrate,\
               timestep,totaltime,boundary_conditions,BC_velo_dof,BC_velo_value,\
               x_max,x_min,y_max,y_min,z_max,z_min,boundaries,nsvars_beam,nsvars_conn_t,\
               nsvars_conn_l,nsecgp,nsvars_secgp,cellwallthickness_early,merge_operation,merge_tol,\
               precrackFlag,precrack_elem)
    
[mass,volume,density,porosity] = \
               ModelInfo_precrack(box_shape,boundary_points,z_min,z_max,\
               nbeam_per_grain,fiberlength,NURBS_degree,nctrlpt_per_beam,\
               props_beam,props_connector_t_bot,props_connector_t_reg,\
               props_connector_t_top,props_connector_l,iprops_beam,\
               beam_connectivity,connector_t_bot_connectivity,\
               connector_t_reg_connectivity,connector_t_top_connectivity,\
               connector_l_connectivity,ConnMeshData,x_indent_size,y_indent_size)
    
ax.set_title('mass = {:8.4e} kg, bulk density = {:8.4e} kg/m3, porosity = {:8.4e}'.format(mass*1e3,density*1e12,porosity))

plt.savefig(Path('meshes/' + geoName + '/' + geoName + '.png'))

FileTime = time.time() 
print('Files generated in {:.3f} seconds'.format(FileTime - BeamTime))

# ==================================================================
# Generate log file for the mesh generation
LogFile(geoName,iter_max,r_min,r_max,nrings,width_heart,width_early,width_late,\
        log_center,box_shape,box_center,box_size,x_min,x_max,y_min,y_max,
        cellsize_early,cellsize_late,cellwallthickness_early,cellwallthickness_late,\
        merge_operation,merge_tol,precrackFlag,precrack_widths,boundaryFlag,\
        fiberlength,theta_min,theta_max,z_min,z_max,long_connector_ratio,\
        NURBS_degree,nctrlpt_per_beam,nconnector_t_precrack,nconnector_l_precrack,\
        nParticles,nbeamElem,mass,volume,density,porosity,\
        stlFlag,inpFlag,inpType,radial_growth_rule,\
        startTime,placementTime,voronoiTime,RebuildvorTime,BeamTime,FileTime)