import math as math
import numpy as np
from scipy.spatial.distance import cdist

def rotate_around_point_highperf(xy, radians, origin=(0, 0)):
    """Rotate a point around a given point.
    
    I call this the "high performance" version since we're caching some
    values that are needed >1 time. It's less readable than the previous
    function but it's faster.
    """
    x, y = xy
    offset_x, offset_y = origin
    adjusted_x = (x - offset_x)
    adjusted_y = (y - offset_y)
    cos_rad = math.cos(radians)
    sin_rad = math.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y

    return qx, qy

NURBS_degree = 2
z_min = 0
z_max = 20.0
grain_length = 1 # mm
height = z_max-z_min

# nsegments = 10
nsegments = int(height/grain_length)

theta_min = 0.0
theta_max = 0.0

long_connector_ratio = 0.2
npt_per_layer = 131
generation_center = (0,0)
max_wings = 4
nridge = 196
nvertex = 131

voronoi_vertices = np.array([[-1.00000000e+01, -1.00000000e+01],
       [-1.00000000e+01, -2.61368802e+00],
       [-1.00000000e+01,  2.68371152e+00],
       [-1.00000000e+01,  3.84147211e+00],
       [-1.00000000e+01,  6.50952213e+00],
       [-1.00000000e+01,  8.98308519e+00],
       [-1.00000000e+01,  1.00000000e+01],
       [-9.68439793e+00,  1.00000000e+01],
       [-9.57859391e+00,  9.47363731e+00],
       [-9.41876226e+00,  3.60330134e+00],
       [-8.85769382e+00, -3.14586129e+00],
       [-8.83266883e+00, -1.00000000e+01],
       [-8.70394761e+00, -3.38239181e+00],
       [-8.64259083e+00, -6.89031612e+00],
       [-8.36852334e+00,  4.51315344e+00],
       [-8.36041798e+00,  4.58379417e+00],
       [-8.03825511e+00, -1.00000000e+01],
       [-7.85650180e+00,  1.00000000e+01],
       [-7.47261301e+00,  9.25424661e+00],
       [-7.17552769e+00, -1.68704001e+00],
       [-6.62474609e+00,  3.18927576e+00],
       [-6.50320822e+00,  8.77147880e+00],
       [-6.16365363e+00,  1.00000000e+01],
       [-6.07323787e+00,  8.99631949e-02],
       [-6.06339403e+00, -9.41457542e+00],
       [-6.04275470e+00,  7.47373130e+00],
       [-5.47687786e+00,  7.46267237e+00],
       [-5.46715848e+00, -1.00000000e+01],
       [-4.82765232e+00,  3.05665738e+00],
       [-4.16879663e+00, -8.21073602e+00],
       [-3.88117689e+00, -5.07218200e+00],
       [-3.69197987e+00,  4.45818618e-02],
       [-3.68506477e+00,  6.54686311e-02],
       [-3.49737509e+00,  3.41441979e+00],
       [-3.48707983e+00, -6.32040315e+00],
       [-3.37459541e+00, -6.60665951e+00],
       [-2.73603342e+00,  4.55315329e+00],
       [-2.66647890e+00, -1.90788253e+00],
       [-2.64799833e+00, -1.00000000e+01],
       [-2.24062286e+00,  1.72764872e+00],
       [-2.20758714e+00, -2.35173791e+00],
       [-1.73971402e+00, -2.26311712e+00],
       [-1.67473033e+00,  4.93581745e+00],
       [-1.35287157e+00, -6.92881144e+00],
       [-1.34127760e+00, -2.32659161e+00],
       [-1.25990865e+00,  1.59584050e+00],
       [-1.20119724e+00,  4.32253432e-01],
       [-1.04067775e+00,  1.00000000e+01],
       [-7.53721214e-01,  1.88491373e+00],
       [-5.50305544e-01, -2.50231271e+00],
       [-4.57678236e-01, -1.00000000e+01],
       [-3.63679172e-01,  2.16319555e-01],
       [-3.11968565e-01,  9.50566123e+00],
       [-2.94093064e-01,  1.92853832e+00],
       [-2.91394263e-01, -2.72121456e+00],
       [-2.17553648e-01, -1.80011236e-02],
       [-1.61148895e-01,  3.00606606e-01],
       [-3.94955966e-02,  2.15353290e-01],
       [-2.06874988e-03,  9.88592321e-03],
       [ 7.24214032e-03,  3.26810936e-02],
       [ 1.81146354e-02,  2.51687771e-01],
       [ 9.00000000e-02,  1.11000000e+00],
       [ 1.47191701e-01,  3.99713111e-01],
       [ 1.56538703e-01, -1.60533689e+00],
       [ 2.49568194e-01,  1.46484155e-01],
       [ 2.51826786e-01,  1.68836530e-01],
       [ 2.65745325e-01,  2.18577395e-01],
       [ 4.83298608e-01, -1.00998356e+00],
       [ 5.42561853e-01, -6.24604836e+00],
       [ 7.61396242e-01,  4.83288528e-01],
       [ 8.07030645e-01, -5.93827187e-01],
       [ 1.15904452e+00, -6.43228216e+00],
       [ 1.24076094e+00,  2.19321409e-01],
       [ 1.38697600e+00, -5.78328826e-02],
       [ 1.47183379e+00,  2.93299476e+00],
       [ 1.66546433e+00, -6.19542017e+00],
       [ 1.77614405e+00,  9.30330397e+00],
       [ 1.79108312e+00,  7.77191508e+00],
       [ 2.07295450e+00, -2.91426397e+00],
       [ 2.10915596e+00, -1.00000000e+01],
       [ 2.16161747e+00, -2.89051564e+00],
       [ 2.29169266e+00,  3.67978867e+00],
       [ 2.30257256e+00,  1.00000000e+01],
       [ 2.30563712e+00, -3.01802011e+00],
       [ 2.75548454e+00, -5.98629183e+00],
       [ 2.78065935e+00,  5.44377647e+00],
       [ 3.08229931e+00, -1.60039700e+00],
       [ 3.16176709e+00, -9.34364992e-01],
       [ 3.23501177e+00, -2.73557975e-01],
       [ 3.30665974e+00, -4.08204321e-01],
       [ 3.30803296e+00,  1.08067251e+00],
       [ 3.33733926e+00, -1.80935681e+00],
       [ 3.55245777e+00, -1.00000000e+01],
       [ 3.65507101e+00, -2.20225847e+00],
       [ 3.91825433e+00, -5.84349995e+00],
       [ 4.09302398e+00, -9.56889666e+00],
       [ 4.38782605e+00, -5.50386805e+00],
       [ 4.46363622e+00, -3.94900173e+00],
       [ 4.48045351e+00, -4.63156453e+00],
       [ 4.63536795e+00, -4.92366719e+00],
       [ 4.65685633e+00,  1.33418496e+00],
       [ 4.76771287e+00, -2.69182861e+00],
       [ 5.05731891e+00, -2.25541911e+00],
       [ 5.39953969e+00, -4.61554142e+00],
       [ 5.53887533e+00,  4.22209810e+00],
       [ 6.00013294e+00, -1.57329407e+00],
       [ 6.36583628e+00,  3.22439102e+00],
       [ 6.55500000e+00, -8.10000000e-01],
       [ 6.62278055e+00,  2.80210738e+00],
       [ 6.87944363e+00, -7.42673342e-01],
       [ 6.94281681e+00,  1.00000000e+01],
       [ 7.11645971e+00, -3.16585890e-02],
       [ 7.44264671e+00,  2.01252432e+00],
       [ 7.45375068e+00,  1.86075075e+00],
       [ 7.86762706e+00, -8.75693414e+00],
       [ 8.47771846e+00,  7.55340076e+00],
       [ 8.66353285e+00,  6.72671223e+00],
       [ 8.69641186e+00,  6.69190367e+00],
       [ 8.93698892e+00, -1.00000000e+01],
       [ 9.26500847e+00, -8.96620860e+00],
       [ 9.32603582e+00,  5.37352356e+00],
       [ 1.00000000e+01, -1.00000000e+01],
       [ 1.00000000e+01, -8.47760627e+00],
       [ 1.00000000e+01, -4.14089276e+00],
       [ 1.00000000e+01,  9.37515436e-02],
       [ 1.00000000e+01,  6.54279794e-01],
       [ 1.00000000e+01,  3.33153296e+00],
       [ 1.00000000e+01,  4.48749960e+00],
       [ 1.00000000e+01,  6.70677167e+00],
       [ 1.00000000e+01,  7.94462028e+00],
       [ 1.00000000e+01,  1.00000000e+01]])


voronoi_ridges = np.array([[ 30,  34],
       [ 30,  12],
       [ 34,  13],
       [ 13,  12],
       [ 32,  39],
       [ 32,  28],
       [ 39,  33],
       [ 28,  33],
       [ 37,  31],
       [ 37,  19],
       [ 31,  23],
       [ 19,  23],
       [ 32,  31],
       [ 46,  41],
       [ 46,  45],
       [ 40,  41],
       [ 40,  37],
       [ 39,  45],
       [ 40,  30],
       [ 10,  12],
       [ 10,  19],
       [ 20,  28],
       [ 20,  23],
       [  9,  14],
       [ 20,  14],
       [ 74,  53],
       [ 74,  72],
       [ 53,  61],
       [ 61,  69],
       [ 72,  69],
       [ 61,  62],
       [ 62,  66],
       [ 69,  66],
       [ 54,  78],
       [ 54,  49],
       [ 49,  63],
       [ 63,  87],
       [ 87,  86],
       [ 86,  80],
       [ 80,  78],
       [ 75,  84],
       [ 75,  78],
       [ 83,  84],
       [ 83,  80],
       [107, 100],
       [107,  89],
       [ 90,  88],
       [ 90, 100],
       [ 88,  89],
       [ 81,  74],
       [ 81,  90],
       [ 73,  72],
       [ 73,  88],
       [ 67,  70],
       [ 67,  63],
       [ 73,  70],
       [ 87,  89],
       [101,  97],
       [101, 102],
       [102, 105],
       [105, 103],
       [103,  99],
       [ 97,  98],
       [ 99,  98],
       [ 93,  97],
       [ 93, 101],
       [ 91, 102],
       [ 91,  93],
       [ 91,  86],
       [107, 105],
       [123, 109],
       [123, 103],
       [109, 107],
       [ 83,  98],
       [ 25,  15],
       [ 25,  21],
       [ 21,  18],
       [ 18,   8],
       [ 14,  15],
       [ 26,  36],
       [ 26,  52],
       [ 52,  76],
       [ 36,  42],
       [ 42,  77],
       [ 76,  77],
       [ 33,  36],
       [ 25,  26],
       [ 45,  48],
       [ 48,  42],
       [ 81,  85],
       [ 48,  53],
       [ 85,  77],
       [ 44,  41],
       [ 44,  55],
       [ 46,  51],
       [ 51,  55],
       [ 61,  56],
       [ 51,  56],
       [ 44,  49],
       [ 67,  58],
       [ 55,  58],
       [111, 113],
       [109, 111],
       [108, 112],
       [108, 100],
       [112, 113],
       [114,  95],
       [114,  96],
       [ 95,  94],
       [ 96,  94],
       [119, 114],
       [ 99,  96],
       [ 84,  94],
       [106, 120],
       [106, 104],
       [104, 116],
       [120, 117],
       [117, 116],
       [108, 106],
       [ 85, 104],
       [116, 115],
       [ 59,  57],
       [ 59,  64],
       [ 65,  64],
       [ 65,  60],
       [ 60,  57],
       [ 56,  57],
       [ 58,  59],
       [ 70,  64],
       [ 66,  65],
       [ 62,  60],
       [ 43,  68],
       [ 68,  71],
       [ 75,  71],
       [ 34,  35],
       [ 54,  68],
       [ 35,  43],
       [ 24,  29],
       [ 35,  29],
       [  9,   3],
       [  9,   2],
       [  1,  10],
       [ 11,  13],
       [122, 119],
       [119, 118],
       [ 15,   4],
       [  8,   5],
       [ 17,  18],
       [  7,   8],
       [ 21,  22],
       [ 47,  52],
       [ 82,  76],
       [111, 124],
       [125, 113],
       [112, 126],
       [127, 120],
       [128, 117],
       [110, 115],
       [129, 115],
       [ 50,  43],
       [ 79,  71],
       [ 95,  92],
       [ 38,  29],
       [ 27,  24],
       [ 16,  24],
       [  3,   2],
       [  2,   4],
       [  4,   5],
       [  5,   6],
       [  6,   7],
       [  7,  17],
       [ 17,  22],
       [ 22,  47],
       [ 47,  82],
       [ 82, 110],
       [110, 130],
       [130, 129],
       [129, 128],
       [128, 127],
       [127, 126],
       [126, 125],
       [125, 124],
       [124, 123],
       [123, 122],
       [122, 121],
       [121, 118],
       [118,  92],
       [ 92,  79],
       [ 79,  50],
       [ 50,  38],
       [ 38,  27],
       [ 27,  16],
       [ 16,  11],
       [ 11,   0],
       [  0,   1],
       [  1,   3]], 'int64')

flattened_all_vertices_2D = np.array([[ 1.93000000e+02,  1.94000000e+02,  1.41000000e+02, ...,
         1.77000000e+02,  1.75000000e+02,  1.76000000e+02],
       [ 1.10000000e+01,  1.00000000e+00,  1.00000000e+01, ...,
         1.28000000e+02,  1.10000000e+02,  1.29000000e+02],
       [ 5.83665585e-01,  3.69315599e+00,  6.30093605e-01, ...,
         6.18924308e-01,  1.52859159e+00,  1.02768986e+00],
       [ 6.00000000e-03,  6.00000000e-03,  6.00000000e-03, ...,
         6.00000000e-03,  6.00000000e-03,  6.00000000e-03],
       [ 6.28318531e+00,  1.57079633e+00, -4.35977859e-01, ...,
        -1.57079633e+00,  3.14159265e+00, -1.57079633e+00],
       [-1.00000000e+00,  1.00000000e+00,  1.00000000e+00, ...,
         1.00000000e+00, -1.00000000e+00,  1.00000000e+00]])



nctrlpt_per_elem = NURBS_degree + 1
nctrlpt_per_beam = 2*NURBS_degree + 1

nlayers = nctrlpt_per_beam*nsegments

nconnector_t_per_beam = int((nctrlpt_per_beam-1)/NURBS_degree+1)
nconnector_t_per_grain = int(nconnector_t_per_beam*nsegments)

segment_length = (z_max - z_min) / (nsegments + (nsegments-1)*long_connector_ratio)
segment_angle = (theta_max-theta_min) / (nsegments + (nsegments-1)*long_connector_ratio)

connector_l_angle = segment_angle*long_connector_ratio
connector_l_length = segment_length*long_connector_ratio

# theta = np.linspace(theta_min,theta_max,nlayers-(nsegments-1))
# z_coord = np.linspace(z_min,z_max,nlayers-(nsegments-1))

theta = np.linspace(theta_min,theta_max-connector_l_angle*(nsegments-1),nlayers-(nsegments-1))
z_coord = np.linspace(z_min,z_max-connector_l_length*(nsegments-1),nlayers-(nsegments-1))

# Insert repeated layers for the longitudinal connectors
for i in range(nlayers-(nsegments-1)-(nctrlpt_per_beam-1),1,-(nctrlpt_per_beam-1)): 
    theta = np.insert(theta,i,theta[i-1])
    theta[i:] += connector_l_angle
for i in range(nlayers-(nsegments-1)-(nctrlpt_per_beam-1),1,-(nctrlpt_per_beam-1)):
    z_coord = np.insert(z_coord,i,z_coord[i-1])
    z_coord[i:] += connector_l_length

# Data for IGA use
vertices_new = np.zeros((nlayers,npt_per_layer,3))
for i in range(0,nlayers):
    for j in range(0,npt_per_layer):
        vertices_new[i,j,:2] = rotate_around_point_highperf(voronoi_vertices[j,:], theta[i], generation_center)
        ### add noise here
        vertices_new[i,j,0] = vertices_new[i,j,0]*0.01*np.random.random()
        vertices_new[i,j,1] = vertices_new[i,j,1]*0.01*np.random.random()
        vertices_new[i,j,2] = z_coord[i]
        

# Vertex Data for IGA
IGAvertices = np.reshape(vertices_new,(-1,3))

# Connectivity for IGA
npt_total = nlayers*npt_per_layer
# vertex_connectivity = np.linspace(0,npt_total-1,npt_total)

# Beams
ngrain = nvertex
nbeam_total = ngrain*nsegments
beam_connectivity_original = np.zeros((nbeam_total,nctrlpt_per_beam))
for i in range(0,nsegments):
    for j in range(0, ngrain):
        irow = i*ngrain + j
        ivertex = i*npt_per_layer*nctrlpt_per_beam + j
        for icol in range(0,nctrlpt_per_beam):
            beam_connectivity_original[irow,icol] = ivertex+icol*npt_per_layer

# Rearange beam connectivity such that each row is corresponding to a Bezier beam element
beam_connectivity = np.copy(beam_connectivity_original)
for ictrlpt in range(nctrlpt_per_beam-NURBS_degree,1,-NURBS_degree):
    beam_connectivity = np.insert(beam_connectivity,ictrlpt,beam_connectivity[:,ictrlpt-1],axis=1)

beam_connectivity_original = (beam_connectivity_original+1).astype(int) # +1 because of in abaqus index starts from 1
beam_connectivity = (np.reshape(beam_connectivity,(-1,nctrlpt_per_elem))+1).astype(int) # +1 because of in abaqus index starts from 1


nbeamElem = beam_connectivity.shape[0]


# Transverse connectors 
nconnector_t = nridge*nconnector_t_per_grain
connector_t_connectivity = np.zeros((nconnector_t,2))
for i in range(0,nsegments):
    for j in range(0,nconnector_t_per_beam):
        for k in range(0,nridge):
            irow = i*nconnector_t_per_beam*nridge + j*nridge + k
            connector_t_connectivity[irow,:] = voronoi_ridges[k]+i*npt_per_layer*nctrlpt_per_beam+j*NURBS_degree*npt_per_layer
            
connector_t_connectivity = (connector_t_connectivity+1).astype(int)
# Slicing block indices, https://stackoverflow.com/questions/39692769/efficient-numpy-indexing-take-first-n-rows-of-every-block-of-m-rows
connector_t_index = np.linspace(0,nconnector_t-1,nconnector_t).astype(int)
connector_t_bot_index = connector_t_index[np.mod(np.arange(connector_t_index.size),nridge*nconnector_t_per_beam)<nridge]
connector_t_top_index = connector_t_index[np.mod(np.arange(connector_t_index.size),nridge*nconnector_t_per_beam)>=(nconnector_t_per_beam-1)*nridge]
connector_t_reg_index = np.setdiff1d(np.setdiff1d(connector_t_index, connector_t_bot_index),connector_t_top_index)
connector_t_bot_connectivity = np.copy(connector_t_connectivity)[connector_t_bot_index,:]
connector_t_top_connectivity = np.copy(connector_t_connectivity)[connector_t_top_index,:]
connector_t_reg_connectivity = np.copy(connector_t_connectivity)[connector_t_reg_index,:]


# Longitudinal connectors
# nwings = all_vertices_2D[:,2].astype(int)
nwings = np.array([2,3,3,3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
           3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
             3,3,3,3,3,3,2,3,3,3,3,3,3,2])
nconnector_l_per_layer = sum(nwings)
nconnector_l = nconnector_l_per_layer * (nsegments-1)
connector_l_connectivity = np.zeros((nconnector_l,2))
connector_l_vertex_dict = np.zeros(nconnector_l)
irow_conn_l = 0
for i in range(0,nsegments-1): # loop over layers of longitudinal connectors
    for j in range(0,ngrain): # loop over grains in each layer
        n = nwings[j]
        irow = i*ngrain + j
        for k in range(0,n): # loop over wings of each grain
            connector_l_connectivity[irow_conn_l,:] = (beam_connectivity_original[irow,-1],beam_connectivity_original[irow,-1]+npt_per_layer)
            connector_l_vertex_dict[irow_conn_l] = irow
            irow_conn_l += 1

connector_l_vertex_dict = connector_l_vertex_dict.astype(int)
connector_l_connectivity = connector_l_connectivity.astype(int)

nconnector_total = nconnector_t + nconnector_l
    



    # return IGAvertices,vertex_connectivity,beam_connectivity_original,nbeam_total,\
    #     beam_connectivity,nbeamElem,nctrlpt_per_beam,nlayers,segment_length,connector_t_connectivity,\
    #     connector_t_bot_connectivity,connector_t_top_connectivity,\
    #     connector_t_reg_connectivity,connector_l_connectivity,nconnector_t_per_beam,\
    #     nconnector_t_per_grain,nconnector_t,nconnector_l,nconnector_total,\
    #     theta,z_coord,connector_l_vertex_dict


# all_vertices_2D = np.array([[-10.        , -10.        ,   2.        , ...,   1.57079633,
#           0.        ,   0.        ],
#        [-10.        ,  -2.61368802,   3.        , ...,   4.71238898,
#           1.57079633,   0.        ],
#        [-10.        ,   2.68371152,   3.        , ...,   1.57079633,
#           1.57079633,   0.        ],
#        ...,
#        [ 10.        ,   6.70677167,   3.        , ...,   1.57079633,
#          -1.57079633,   0.        ],
#        [ 10.        ,   7.94462028,   3.        , ...,   1.57079633,
#          -1.57079633,   0.        ],
#        [ 10.        ,  10.        ,   2.        , ...,  -1.57079633,
#           0.        ,   0.        ]])


# all_ridges = np.array([[ 30,  34, 131, 270, 409, 548],
#        [ 30,  12, 132, 271, 410, 549],
#        [ 34,  13, 133, 272, 411, 550],
#        ...,
#        [ 11,   0, 741, 798, 855, 912],
#        [  0,   1, 742, 799, 856, 913],
#        [  1,   3, 743, 800, 857, 914]], 'int64')