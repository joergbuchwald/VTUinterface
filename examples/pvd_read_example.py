'''
This example demonstrates the usage of vtuIO.PVDIO at the results of the OGS-benchmark
Elliptic problem with Dirichlet-type boundary conditions,
i.e. it requires the presence of:

square_1e2_pcs_0.pvd
square_1e2_pcs_0_ts_0_t_0.000000.vtu
square_1e2_pcs_0_ts_1_t_1.000000.vtu

The pressure at a point is read and plotted over time (two time points)
'''

import matplotlib.pyplot as plt 	# for fancy plots

import vtuIO	# to read and process (point interpolation) vtu- and pvd-files 
# class methods for information
#	PVDIO
#    __init__(self, folder, filename, dim=3):
#    readPVD(self,filename):
#    readTimeSeries(self,fieldname, pts = {'pt0': (0.0,0.0,0.0)}):
#    readTimeStep(self, timestep, fieldname):


# read pvd-file specified by path and filename
# dim refers to the actual dimension:
# 2D data in 2D are OK (dim=2), 3D data in 3D are OK (dim=3).
# Note that for 2D data in 3D, e.g. x,y!=0 and z=0, 
# dim must be set to 2, otherwise the interpolator fails.
# Currently PVDIO assumes 2D data in 3D at x,y and ignores z.
pvdfile=vtuIO.PVDIO(".", "square_1e2_pcs_0.pvd", dim=2)

# get time vector from pvd-data (dictionary)
time=pvdfile.ts_files['ts']   

# define points for interpolation (dictionary)
selected_points={'pt0': (0.25, 0.5, 0.0), 'pt1': (0.75, 0.5, 0.0)}

# read and interpolate from vtu-files listed in pvd
pressure_interpolation=pvdfile.readTimeSeries('pressure', selected_points)

# read pressure at pt0 from interpolations (dictionary)
pressure_pt0=pressure_interpolation['pt0']    

# plot some result
plt.plot(time, pressure_pt0)
titlestring="At point "+str(selected_points['pt0'])
plt.title(titlestring)
plt.xlabel('t')
plt.ylabel('p')
plt.show()

# do something with pt1 or whatever you like
# ...

