def inputParams():

    geoName = 'Debug_QML_v11_010824'

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

    return geoName, radial_growth_rule, iter_max, print_interval, \
        r_min, r_max, nrings, width_heart, width_early, width_late, generation_center, \
        cellsize_early, cellsize_late, cellwallthickness_early, cellwallthickness_late, \
        boundaryFlag, box_shape, box_center, box_size, \
        nsegments, theta_max, theta_min, z_max, z_min, long_connector_ratio, \
        skeleton_density, merge_operation, merge_tol, precrack_widths, precrackFlag, \
        stlFlag, inpFlag, inpType