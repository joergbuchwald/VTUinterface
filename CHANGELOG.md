# Changelog

All notable changes to **VTUinterface** will be documented in this file.

## [0.681]

### Changes
* new tool that enables spatial transformation of slices
* more functionalities to read and interpolate cell data based on cell center points
* VTUinterface is now able to read pvtu files

## [0.68]

### Changes
* changes in interface to distinguish between cell data an point data
* filename argument of PVDIO contains now the directory argument as well (the folder keyword argument is dropped)
* more tests
* VTUinterface can deal with different orientations for 1d and 2d
* VTK backend is added to enable better interpolation as voronoi based interpolation in scipy often fails

### Bugfixes

### Additions

