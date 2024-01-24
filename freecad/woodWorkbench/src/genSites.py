import freecad.woodWorkbench.tools.WoodMeshGenTools_v11 as WoodMeshGen
from pathlib import Path
import os

def genSites(dumForm,outdir):

    [geoName, radial_growth_rule, iter_max, print_interval, \
        r_min, r_max, nrings, width_heart, width_early, width_late, generation_center, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_size, \
        nsegments, theta_max, theta_min, z_max, z_min, long_connector_ratio, \
        skeleton_density, merge_operation, merge_tol, precrack_widths, precrackFlag, \
        stlFlag, inpFlag, inpType] = dumForm

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
        


    return sites, radii