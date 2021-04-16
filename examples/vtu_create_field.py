'''
This example demonstrates the usage of vtuIO.VTUIO to create point field (arrays)
on the mesh.
'''


import vtuIO	# to read and process (point interpolation) vtu- and pvd-files 


# read file
data=vtuIO.VTUIO("square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)

def getPressure(x,y,z):
    if x<=0.5:
        return -5.0e3
    else:
        return 5.0e3

data.func2Field(getPressure, "p0","p0_field.vtu")

# multidimensional example

def fct(x,y,z):
    if x<=0.5:
        return -5.0e3*0.3*0.371163
    else:
        return 5.0e3*0.3*0.95

def fct2(x,y,z):
    return 0

# result is a four-dimensional point field
data.func2mdimField([fct,fct,fct2,fct2], "sigma0","sigma0_field.vtu")

