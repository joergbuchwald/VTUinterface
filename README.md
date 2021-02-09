# VTUinterface 

VTUinterface is a python package for easy accessing VTU/PVD files as outputed by Finite Element software like OpenGeoSys. It uses the VTK python wrapper and linear interpolation between time steps and grid points access any points in and and time within the simulation domain.
While beeing a python package, it was also tested in Julia, where it can be accessed via PyCall:


```julia
ENV["PYTHON"] = "/usr/bin/python3"
using Pkg
Pkg.build("PyCall")
```

    [32m[1m   Building[22m[39m Conda ─→ `~/.julia/packages/Conda/x5ml4/deps/build.log`
    [32m[1m   Building[22m[39m PyCall → `~/.julia/packages/PyCall/tqyST/deps/build.log`



```julia
using PyCall
@pyimport vtuIO
```

VTUinterface together with ogs6py can be viewed in action here:

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/eihNKjK-I-s/0.jpg)](https://www.youtube.com/watch?v=eihNKjK-I-s)

Single VTU files can be accessed via:

# 1. reading a single VTU file


```julia
vtufile = vtuIO.VTUIO("examples/square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)
```




    PyObject <vtuIO.VTUIO object at 0x7f08f5dfd310>



The `dim` argument is needed for correct interpolation. By defualt `dim=3` is assumed.
Basic VTU properties, like fieldnames, points and the corresponding as provided by the unstructured grid VTK class: 


```julia
fields=vtufile.getFieldnames()
```




    4-element Array{String,1}:
     "D1_left_bottom_N1_right"
     "Linear_1_to_minus1"
     "pressure"
     "v"




```julia
vtufile.points
```




    121×2 Array{Float64,2}:
     0.0  0.0
     0.1  0.0
     0.2  0.0
     0.3  0.0
     0.4  0.0
     0.5  0.0
     0.6  0.0
     0.7  0.0
     0.8  0.0
     0.9  0.0
     1.0  0.0
     0.0  0.1
     0.1  0.1
     ⋮    
     1.0  0.9
     0.0  1.0
     0.1  1.0
     0.2  1.0
     0.3  1.0
     0.4  1.0
     0.5  1.0
     0.6  1.0
     0.7  1.0
     0.8  1.0
     0.9  1.0
     1.0  1.0




```julia
vtufile.getField("v")
```




    121×2 Array{Float64,2}:
     2.0   0.0
     2.0   1.62548e-16
     2.0  -9.9123e-16
     2.0  -9.39704e-16
     2.0  -4.08897e-16
     2.0   1.36785e-16
     2.0  -3.23637e-16
     2.0  -2.30016e-16
     2.0  -7.69185e-16
     2.0  -2.27994e-15
     2.0   1.53837e-15
     2.0   3.25096e-16
     2.0  -3.62815e-16
     ⋮    
     2.0  -8.88178e-16
     2.0   0.0
     2.0  -2.22045e-16
     2.0   9.9123e-16
     2.0  -1.2648e-15
     2.0   5.48137e-16
     2.0  -3.89112e-17
     2.0  -2.03185e-17
     2.0  -1.02098e-15
     2.0  -5.03586e-16
     2.0  -3.37422e-15
     2.0   8.88178e-16



Aside basic VTU properties, the field data at any given point can be retrieved:


```julia
points = Dict("pt0"=> (0.5,0.5,0.0), "pt1"=> (0.2,0.2,0.0))
```




    Dict{String,Tuple{Float64,Float64,Float64}} with 2 entries:
      "pt1" => (0.2, 0.2, 0.0)
      "pt0" => (0.5, 0.5, 0.0)




```julia
# Python: points={'pt0': (0.5,0.5,0.0), 'pt1': (0.2,0.2,0.0)} 
```


```julia
point_data = vtufile.getPointData("pressure", pts=points)
```




    Dict{Any,Any} with 2 entries:
      "pt1" => 0.6
      "pt0" => 3.41351e-17



# 2. Writing VTU files
some simple methods also exist for adding new fields to an existing VTU file or save it separately:


```julia
size = length(vtufile.getField("pressure"))
```




    121




```julia
p0 = ones(size) * 1e6;
```


```julia
# Python: size = len(vtufile.getField("pressure"))
# p0 = np.ones() *1.0e6
```


```julia
vtufile.writeField(p0, "initialPressure", "mesh_initialpressure.vtu")
```


```julia

```

A new field can also created from a three-argument function for all space-dimensions:


```julia
function p_init(x,y,z)
    if x < 0.5
        return -0.5e6
    else
        return +0.5e6
    end
end
```




    p_init (generic function with 1 method)




```julia
# Python:
# def p_init(x,y,z):
#    if x<0.5:
#        return -0.5e6
#    else:
#        return 0.5e6
```


```julia
vtufile.func2Field(p_init, "p_init", "mesh_initialpressure.vtu")
```

It is also possible to write multidimensional arrays using a function.


```julia
function null(x,y,z)
    return 0.0
end
```




    null (generic function with 1 method)




```julia
vtufile.func2mdimField([p_init,p_init,null,null], "sigma00","mesh_initialpressure.vtu")
```

# Read PVD files:

See examples/pvd*.py for further details.


```julia

```
