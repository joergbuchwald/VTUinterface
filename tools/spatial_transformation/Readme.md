# Spatial Transformation of Mesh Points and Mesh Data
This script reads a vtu-file and transforms its points as well as its point- and cell-data to a new coordinate system.
Due to a frequent use-case, the new coordinate system is defined by a slice, as it is done for the OGS tool ``VerticalSliceFromLayers`` (this defines the forward transformation).
In this use case a 2D slice in 3D is transformed to 2D for running simulations on it. 
The results are then transformed back to 3D for spatial visualization.

If the _z_-coordinate is required to be zero for simulations, then remember the _z_-coordinate originally in 2D (should be one value for all points up to numerical precision).
Run the forward transformation with the parameter ``-z 0.0`` and take the simulation results back to 3D by the reverse trafo (``-r`` option) passing the original 2D _z_-coordinate via ``-z <FLOAT>``.
