#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transforms a mesh (points and vectorial/tensorial data) to new coordinates.
A frequent use case are 2D-slices in 3D to be transformed into x-y-coordinates
(pre-processing), running OGS and transforming back into x-y-z coordinates
(post-processing).

AUTOMATIC DETECTION
    runfile('vtu_trafo.py', args='-i example4preprocessing.vtu -o forward.vtu -s')
    runfile('vtu_trafo.py', args='-i forward.vtu -o backward.vtu -r')
MANUAL SLICE DEFINITION
    runfile('vtu_trafo.py', args='-i example4preprocessing.vtu -o forward.vtu -x 9200 -X 18000 -y 9000 -Y 20000')
    runfile('vtu_trafo.py', args='-i forward.vtu -o backward.vtu -x 9200 -X 18000 -y 9000 -Y 20000 -r')
EXAMPLE FOR TRAFO OF RESULTS
    runfile('vtu_trafo.py', args='-i example4postprocessing.vtu -o backward.vtu -x 9200 -X 18000 -y 9000 -Y 20000 -r')

TODO
    implement for further THM result types, possibly for input data too
        done Darcy velocity (point data)
"""
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse

tested_vtk_version = "9.0.1"
vtutrafo_version = "0.0"
vector_names = {'DarcyVelocity'}   # list names of data that transforms like vectors
coordinate_init_value = 0.0   # to detect if not set or invalid
ORG_X="original_x"
ORG_Y="original_y"
ORG_Z="original_z"
stored_data_found = False

# parsing command line arguments

parser = argparse.ArgumentParser(
    description="Transforms mesh points and results between two orthonormal coordinate systems. Needed to run OpenGeoSys on 2D-slices from a 3D-model.",
    epilog="Tested with VTK "
    + tested_vtk_version
)

required = parser.add_argument_group('required arguments')
required.add_argument(
    "-i",
    "--input",
    default="",
    required=True,
    help='input file'
)
required.add_argument(
    "-o",
    "--output",
    default="",
    required=True,
    help="output file"
)

parser.add_argument(
    "-x",
    "--start-x",
    type=float,
    default=coordinate_init_value,
    help="start x coordinate of slice (default is zero, if none of slice coordinates is given, then automatic detection)"
)
parser.add_argument(
    "-X",
    "--end-x",
    type=float,
    default=coordinate_init_value,
    help="end x coordinate of slice (default is zero, if none of slice coordinates is given, then automatic detection)"
)
parser.add_argument(
    "-y",
    "--start-y",
    type=float,
    default=coordinate_init_value,
    help="start y coordinate of slice (default is zero, if none of slice coordinates is given, then automatic detection)"
)
parser.add_argument(
    "-Y",
    "--end-y",
    type=float,
    default=coordinate_init_value,
    help="end y coordinate of slice (default is zero, if none of slice coordinates is given, then automatic detection)"
)
parser.add_argument(
        "-z",
        "--set-2D-z",
        metavar="2D_Z",
        type=float,
        help="enforces given z-coordinate in 2D, except for reverse trafo with given original points (autodetect). Note that forward trafo ends in 2D and reverse trafo starts in 2D"
    )


mutex = parser.add_mutually_exclusive_group()
mutex.add_argument(
        "-r",
        "--reverse",
        action="store_true",
        help="reverse trafo"
    )

mutex.add_argument(
        "-s",
        "--store",
        action="store_true",
        help="store orginal data (only for direct trafo)"
    )

parser.add_argument('-v', '--version', action='version', version=vtutrafo_version)
args = parser.parse_args()

# check command line arguments
store_flag = args.store
reverse_flag = args.reverse
start_x = args.start_x
end_x = args.end_x
start_y = args.start_y
end_y = args.end_y
inputfile = args.input
outputfile = args.output

auto_detect = (start_x == end_x) and (start_y == end_y)

# Read file
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(inputfile)
reader.Update()
vtk_mesh = reader.GetOutput()
vtk_pdata = vtk_mesh.GetPointData()
vtk_cdata = vtk_mesh.GetCellData()
vtk_points = vtk_mesh.GetPoints()

print("{} Points".format(vtk_mesh.GetNumberOfPoints()))
pdata_array_names = []
for i in range(vtk_pdata.GetNumberOfArrays()):
    pdata_array_names.append(vtk_pdata.GetArrayName(i))
print("Point Data:", pdata_array_names)

print("{} Cells".format(vtk_mesh.GetNumberOfCells()))
cdata_array_names = []
for i in range(vtk_cdata.GetNumberOfArrays()):
    cdata_array_names.append(vtk_cdata.GetArrayName(i))
print("Cell Data: ", cdata_array_names)

points = vtk_to_numpy(vtk_points.GetData())
N_points = len(points)

# store original points
if args.store:
    print("store original point coordinates")
    r = vtk_pdata.AddArray(numpy_to_vtk(points[:, 0]))
    vtk_pdata.GetArray(r).SetName(ORG_X)
    r = vtk_pdata.AddArray(numpy_to_vtk(points[:, 1]))
    vtk_pdata.GetArray(r).SetName(ORG_Y)
    r = vtk_pdata.AddArray(numpy_to_vtk(points[:, 2]))
    vtk_pdata.GetArray(r).SetName(ORG_Z)

if (ORG_X in pdata_array_names) and (ORG_Y in pdata_array_names) and (ORG_Z in pdata_array_names):
    stored_data_found = True
    org_points = np.stack( [np.array(vtk_to_numpy(vtk_pdata.GetArray(ORG_X))),
                            np.array(vtk_to_numpy(vtk_pdata.GetArray(ORG_Y))),
                            np.array(vtk_to_numpy(vtk_pdata.GetArray(ORG_Z)))], axis=-1 )


###   DEFINE BASES   ###
global_ex = np.array([1, 0, 0])   # 3D coordinate system
global_ey = np.array([0, 1, 0])
global_ez = np.array([0, 0, 1])
global_base = np.matrix([global_ex, global_ey, global_ez])# note, base as row_vectors

# construct base from given slice direction or automatically if not given
if auto_detect:
    if reverse_flag:
        if stored_data_found:
            print("reading original points from " + inputfile)
            # vectors from point with minimum of a coordinate to point with maximum
            dP = np.array([ org_points[np.argmax(org_points[:,0])] - org_points[np.argmin(org_points[:,0])],
                            org_points[np.argmax(org_points[:,1])] - org_points[np.argmin(org_points[:,1])],
                            org_points[np.argmax(org_points[:,2])] - org_points[np.argmin(org_points[:,2])] ])
        else:
            print("WARNING! Automatic detection for reverse trafo, but no stored data found. No output written.")
            sys.exit()
    else:
        print("detecting slice orientation automatically")
    # vectors from point with minimum of a coordinate to point with maximum
        dP = np.array([ points[np.argmax(points[:,0])] - points[np.argmin(points[:,0])],
                        points[np.argmax(points[:,1])] - points[np.argmin(points[:,1])],
                        points[np.argmax(points[:,2])] - points[np.argmin(points[:,2])] ])
    dPnorm = np.linalg.norm(dP, axis=-1)
    dP_sorted = dP[np.argsort(dPnorm)]
    local_Ex = dP_sorted[-1]   # 2D coordinate system
    local_EY = dP_sorted[-2]
else:
    local_Ex = np.array([end_x-start_x, end_y-start_y, 0.0])   # 2D coordinate system
    local_EY = np.array([0.0, 0.0, 1.0])

local_ex = local_Ex / np.linalg.norm(local_Ex)    # unit vector
local_Ey = local_EY - np.dot(local_EY, local_ex)*local_ex   # make orthogonal
local_ey = local_Ey / np.linalg.norm(local_Ey)   # unit vector
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

# transform mesh points (nodes)
if reverse_flag and auto_detect:
    points = org_points
else:
    if (args.set_2D_z is not None) and reverse_flag:
        points[:,2] = args.set_2D_z   # set z coordinate in 2D, i.e. before trafo
    points = spR.apply(points, inverse=reverse_flag)
    if (args.set_2D_z is not None) and not reverse_flag:
        points[:,2] = args.set_2D_z   # set z coordinate in 2D, i.e. after trafo
vtk_points.SetData(numpy_to_vtk(points))

if stored_data_found:
    vtk_pdata.RemoveArray(ORG_X)
    vtk_pdata.RemoveArray(ORG_Y)
    vtk_pdata.RemoveArray(ORG_Z)

# transform point data
for pdata_array_name in pdata_array_names:
    if pdata_array_name in vector_names:
        vtk_pdata_array = vtk_pdata.GetArray(pdata_array_name)
        pdata_array = vtk_to_numpy(vtk_pdata_array)
        rows, cols = pdata_array.shape
        if cols == 2:   # append zero z-component
            pdata_array = spR.apply( np.concatenate((pdata_array, np.zeros((rows,1))), axis=1),
                    inverse=reverse_flag)
            print(pdata_array_name + " augmented (2D to 3D) and transformed")
            vtk_pdata.RemoveArray(pdata_array_name)
            r = vtk_pdata.AddArray(numpy_to_vtk(pdata_array))
            vtk_pdata.GetArray(r).SetName(pdata_array_name)
        elif cols == 3:
            pdata_array = spR.apply(pdata_array, inverse=reverse_flag)
            print(pdata_array_name + " transformed")
            vtk_pdata.RemoveArray(pdata_array_name)
            r = vtk_pdata.AddArray(numpy_to_vtk(pdata_array))
            vtk_pdata.GetArray(r).SetName(pdata_array_name)
        else:
            print("Not a vector: ", pdata_array_name)

# write file
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(outputfile)
writer.SetInputData(vtk_mesh)
writer.Write()
