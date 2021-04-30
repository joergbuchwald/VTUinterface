'''
This example demonstrates the usage of vtuIO.VTUIO at the results of the OGS-benchmark
Elliptic problem with Dirichlet-type boundary conditions,
i.e. it requires the presence of:

square_1e2_pcs_0_ts_1_t_1.000000.vtu

Data field names and points are read and printed, 
and finally the data field "pressure" read and printed
'''

import matplotlib.pyplot as plt 	# for fancy plots
import matplotlib.tri as tri    # for triangulation

import vtuIO	# to read and process (point interpolation) vtu- and pvd-files 
# class methods for information
#	VTUIO		
#    __init__(self, filename, dim=3):
#    getNeighbors(self, points_interpol):
#    getData(self, neighbors, points_interpol, fieldname):
#    getField(self, fieldname):
#    getFieldnames(self):
#    writeField(self, field, fieldname, ofilename):


# read file
data=vtuIO.VTUIO("square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)

# see what is inside the file, which fields and points
fields=data.getFieldnames()
print("fields:")
print(fields)

points=data.points
print("points:")
print(points)

pressure_field = data.getField("pressure")
print("pressure at points")
print(pressure_field)


point_data = data.getPointData("pressure", pts={'pt0':(0.5,0.5,0.4)})
print(point_data)

# contour plot

triang=tri.Triangulation(points[:,0],points[:,1])
fig, ax = plt.subplots(ncols=1,figsize=(20,8))
contour = ax.tricontourf(triang, pressure_field)
fig.colorbar(contour,ax=ax,label='p (Pa)')
plt.show()
