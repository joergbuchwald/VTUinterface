[![DOI](https://zenodo.org/badge/282728412.svg)](https://zenodo.org/badge/latestdoi/282728412)


# VTUinterface 

VTUinterface is a python package for easy accessing VTU/PVD files as outputed by Finite Element software like OpenGeoSys. It uses the VTK python wrapper and linear interpolation between time steps and grid points access any points in and and time within the simulation domain.

[Documentation reference](https://joergbuchwald.github.io/VTUinterface-doc)


VTUinterface together with ogs6py can be viewed in action here:

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/eihNKjK-I-s/0.jpg)](https://www.youtube.com/watch?v=eihNKjK-I-s)

# 0. Installation

clone the repository and use pip to install the package

```shell
# git clone https://github.com/joergbuchwald/VTUinterface.git
# cd VTUinterface
# pip install --user .
```

# 1. Quick start

[Basic Usage (python)](https://github.com/joergbuchwald/VTUinterface/blob/master/README_python.md)

Although, a python package, VTUinterface is tested to work through PyCall under julia as well:

[Basic Usage (julia)](https://github.com/joergbuchwald/VTUinterface/blob/master/README_julia.md)


# FAQ/Troubleshooting

# Troubleshooting


As the input data is triangulated with QHull for the linear interpolation it might fail at boundaries or if a wrong input dimension is given.
Possible solutions:

- In order for interpolation to work correctly providing the correct dimension (set via `dim` keyword) of the problem is crucial.
- For some meshes it might help to adjust the number of points taken into account by the triangulation, which can be done using the `nneighbors` keyword. Default value is 20.
- Especially along boundaries, linear interpolation with the QHULL method often fails, this can be resolved bei using nearest neighbor interpolation.


