'''
This example demonstratesthe extraction of point set data using
vtuIO.PVDIO results of the OGS-benchmark
Elliptic problem with Dirichlet-type boundary conditions,
i.e. it requires the presence of:

square_1e2_pcs_0.pvd
square_1e2_pcs_0_ts_0_t_0.000000.vtu
square_1e2_pcs_0_ts_1_t_1.000000.vtu

The pressure along two axes for two different times is read and plotted over r
'''
import numpy as np

import matplotlib.pyplot as plt 	# for fancy plots

import vtuIO	# to read and process (point interpolation) vtu- and pvd-files 


# read pvd-file specified by path and filename
# dim refers to the actual dimension:
# 2D data in 2D are OK (dim=2), 3D data in 3D are OK (dim=3).
# Note that for 2D data in 3D, e.g. x,y!=0 and z=0, 
# dim must be set to 2, otherwise the interpolator fails.
# Currently PVDIO assumes 2D data in 3D at x,y and ignores z.
pvdfile=vtuIO.PVDIO("square_1e2_pcs_0.pvd", dim=2)

# define xaxis and diagonal (list)
xaxis =  [(i,0,0) for i in np.linspace(start=0.0, stop=1.0, num=100)]
diagonal = [(i,i,0) for i in np.linspace(start=0.0, stop=1.0, num=100)]

# define timestep
t1 = 0.2543
t2 = 0.9

# read and interpolate from vtu-files listed in pvd
pressure_xaxis_t1 = pvdfile.read_set_data(t1, 'pressure', data_type="point", pointsetarray=xaxis)
pressure_diagonal_t1 = pvdfile.read_set_data(t1, 'pressure', data_type="point", pointsetarray=diagonal)
pressure_xaxis_t2 = pvdfile.read_set_data(t2, 'pressure', data_type="point", pointsetarray=xaxis)
pressure_diagonal_t2 = pvdfile.read_set_data(t2, 'pressure', data_type="point", pointsetarray=diagonal)

# convert lists to array:
r_x = np.array(xaxis)[:,0]
r_diag = np.sqrt(np.array(diagonal)[:,0]**2+np.array(diagonal)[:,1]**2)


# plot some result
plt.plot(r_x, pressure_xaxis_t1, label='p_x t=t1')
plt.plot(r_diag, pressure_diagonal_t1, label='p_diag t=t1')
plt.plot(r_x, pressure_xaxis_t2, label='p_x t=t2')
plt.plot(r_diag, pressure_diagonal_t2, label='p_diag t=t2')
titlestring="Pressure along x and diagonal"
plt.title(titlestring)
plt.xlabel('r')
plt.ylabel('p')
plt.legend()
plt.show()

# do something with pt1 or whatever you like
# ...

