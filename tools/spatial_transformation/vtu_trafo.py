#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transforms a mesh (points and vectorial/tensorial data) to new coordinates.
A frequent use case are 2D-slices in 3D to be transformed into x-y-coordinates
(pre-processing), running OGS and transforming back into x-y-z coordinates
(post-processing).

TODO
    implement for THM result types: done Darcy velocity (point data)
    
SPYDER
   runfile('vtu_trafo.py', args='-i example4preprocessing.vtu -o forward.vtu -x 9200 -X 18000 -y 9000 -Y 20000')
   runfile('vtu_trafo.py', args='-i example4postprocessing.vtu -o backward.vtu -x 9200 -X 18000 -y 9000 -Y 20000 -r')
"""
import meshio
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse

tested_meshio_version = "4.4.0"
vtutrafo_version = "0.0"
vector_data = {'DarcyVelocity'}   # list data that transforms like vectors

# parsing command line arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(
    description="Transforms mesh points and results between two orthonormal coordinate systems. Needed to run ogs on 2D-slices from a 3D-model.",
    epilog="Tested with Meshio "
    + tested_meshio_version
)
parser.add_argument(
    "-i",
    "--input",
    default="",
    help='input file'
)
parser.add_argument(
    "-o",
    "--output",
    default="",
    help="output file"
)
parser.add_argument(
    "-x",
    "--start-x",
    type=float,
    default=0.0,
    help="start x coordinate of slice"
)
parser.add_argument(
    "-X",
    "--end-x",
    type=float,
    default=0.0,
    help="end x coordinate of slice"
)
parser.add_argument(
    "-y",
    "--start-y",
    type=float,
    default=0.0,
    help="start y coordinate of slice"
)
parser.add_argument(
    "-Y",
    "--end-y",
    type=float,
    default=0.0,
    help="end y coordinate of slice"
)
parser.add_argument(
        "-r",
        "--reverse",
        action="store_true",
        help="reverse trafo"
    )
parser.add_argument(
        "-z",
        "--set-2D-z",
        type=float,
        help="enforces given z-coordinate in 2D, note that forward trafo ends in 2D and reverse trafo starts in 2D"
    )


parser.add_argument('-v', '--version', action='version', version=vtutrafo_version) 
args = parser.parse_args()

# check command line arguments 
reverse_flag = args.reverse
start_x = args.start_x
end_x = args.end_x
start_y = args.start_y
end_y = args.end_y
inputfile = args.input
outputfile = args.output

###   DEFINE BASES   ###
global_ex = np.array([1, 0, 0])   # 3D coordinate system
global_ey = np.array([0, 1, 0])
global_ez = np.array([0, 0, 1])
global_base = np.matrix([global_ex, global_ey, global_ez])# note, base as row_vectors 

# construct base from slice direction
local_Ex = np.array([end_x-start_x, end_y-start_y, 0.0])   # 2D coordinate system
local_ex = local_Ex / np.linalg.norm(local_Ex)
local_EY = np.array([0.0, 0.0, 1.0])   
local_Ey = local_EY - np.dot(local_EY, local_ex)*local_ex   # make orthogonal
local_ey = local_Ey / np.linalg.norm(local_Ey)
local_ez = np.cross(local_ex, local_ey)
local_base = np.matrix([local_ex, local_ey, local_ez])


###   GENERATE TRANSFORMATION   ###
# my way (to verify scipy way)
#myRmatrix = global_base*local_base.I
#myR = R.from_matrix(myRmatrix)

# scipy way
spR, spMeanError = R.align_vectors(global_base, local_base)   
#spRmatrix = spR.as_matrix()


###   APPLY TRANSFORMATION   ###
mesh=meshio.read(inputfile)   
print(mesh)
points, point_data = mesh.points, mesh.point_data

# transform mesh points (nodes)
if (args.set_2D_z is not None) and reverse_flag:
    points[:,2] = args.set_2D_z   # set z coordinate in 2D, i.e. before trafo
    
new_points = spR.apply(points, inverse=reverse_flag) 

if (args.set_2D_z is not None) and not reverse_flag:
    new_points[:,2] = args.set_2D_z   # set z coordinate in 2D, i.e. after trafo 

mesh.points = new_points 

for point_data_key in point_data.keys():
    if point_data_key in vector_data:
        rows, cols = point_data[point_data_key].shape
        if cols == 2:
            vectors = np.concatenate((point_data[point_data_key], np.zeros((rows,1))), axis=1)
            mesh.point_data[point_data_key] = spR.apply(vectors, inverse=reverse_flag) 
        elif cols == 3:
            vectors = point_data[point_data_key] 
            mesh.point_data[point_data_key] = spR.apply(vectors, inverse=reverse_flag) 
        else:
            print("Not a vector: ", point_data_key)
        
meshio.write(outputfile, mesh)

#pressure_field = data.get_field("pressure")
#print("pressure at points")
#print(pressure_field)
