import freecad.woodWorkbench.tools.WoodMeshGenTools_v11 as WoodMeshGen
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

def clipBox(box_shape,box_center,box_size,boundaryFlag,sites):
    x_min,x_max,y_min,y_max,boundaries,boundary_points,boundarylines = \
        WoodMeshGen.Clipping_Box(box_shape,box_center,box_size,boundaryFlag)
    # *** boundary lines overwritting?

    # if precracked - make inputs eventually ***
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



    return x_min, x_max, y_min, y_max, \
        x_indent, y_indent_min, y_indent_max, x_precrack, y_precrack, \
        vor, boundaries, boundarylines, boundary_points