# -*- coding: utf-8 -*-
"""
VTUinterface is a python package for easy accessing VTU/PVD files as
outputed by Finite Element software like OpenGeoSys. It uses the VTK python
wrapper and linear interpolation between time steps and grid points access
any points in and and time within the simulation domain.

Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
            Distributed under a Modified BSD License.
              See accompanying file LICENSE or
              http://www.opengeosys.org/project/license

"""

# pylint: disable=C0103, R0902, R0914, R0913
import os
import warnings

import numpy as np
import pandas as pd


warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk

from lxml import etree as ET

from scipy.interpolate import griddata
from scipy.interpolate import interp1d


class VTUIO:
    """
    Class for handling I/O of VTU files

    Parameters
    ----------
    filename : `str`
    nneighbors : `int`, optional
                 default: 20
    dim : `int`, optional
          default: 3
    one_d_axis : `int`
                 between 0 and 2, default: 0
    two_d_planenormal : `int`
                 between 0 and 2, default: 2
    interpolation_backend : `str`
                 scipy or vtk
    """
    def __init__(self, filename, nneighbors=20, dim=3, one_d_axis=0, two_d_planenormal=2,
                                                        interpolation_backend="scipy"):
        self.filename = filename
        if filename.split(".")[-1] == "vtu":
            self.reader = vtk.vtkXMLUnstructuredGridReader()
        elif filename.split(".")[-1] == "pvtu":
            self.reader = vtk.vtkXMLPUnstructuredGridReader()
        else:
            raise RuntimeError(f"Did not recognize file extension")
        if os.path.isfile(self.filename) is True:
            self.reader.SetFileName(self.filename)
        else:
            raise RuntimeError(f"File not found: {self.filename}")
        self.reader.Update()
        try:
            self.output = self.reader.GetOutput()
            self.pdata = self.output.GetPointData()
            self.cdata = self.output.GetCellData()
        except AttributeError:
            print(f"File {self.filename} does not contain any data")
        try:
            self.points = vtk_to_numpy(self.output.GetPoints().GetData())
        except AttributeError:
            print(f"File {self.filename} does not contain any points")
        self._cell_center_points = None
        self.dim = dim
        self.nneighbors = nneighbors
        self.one_d_axis=one_d_axis
        self.two_d_planenormal = two_d_planenormal
        if self.dim == 1:
            self.one_d_axis = one_d_axis
            self.points = self.points[:,one_d_axis]
        if self.dim == 2:
            self.plane = [0, 1, 2]
            self.plane.pop(two_d_planenormal)
            self.points = np.delete(self.points, two_d_planenormal, 1)
        self.interpolation_backend = interpolation_backend

    @property
    def header(self):
        header_dict = {"N Cells": [str(len(self.cell_center_points))], "N Points": [str(len(self.points))],
                "N Cell Arrays": [len(self.get_cell_field_names())],
                "N Point Arrays": [len(self.get_point_field_names())],
                "X Bounds": [str(np.min(self.points[:,0])) + ", " + str(np.max(self.points[:,0]))]}
        if self.dim > 1:
            header_dict["Y Bounds"] = [str(np.min(self.points[:,1]))+" "+str(np.max(self.points[:,1]))]
        if self.dim > 2:
            header_dict["Z Bounds"] = [str(np.min(self.points[:,2]))+" "+str(np.max(self.points[:,2]))]
        df = pd.DataFrame(header_dict).T
        return df.rename({0:"Information"}, axis='columns')

    @property
    def data_arrays(self):
        pf_names = self.get_point_field_names()
        cf_names = self.get_cell_field_names()
        data_array_dict = {}
        for name in pf_names:
            field = self.get_point_field(name)
            if field.ndim == 1:
                components = 1
            else:
                components = field.shape[1]
            data_array_dict[name] = ["point", components, np.min(field), np.mean(field), np.max(field)]
        for name in cf_names:
            field = self.get_cell_field(name)
            data_array_dict[name] = ["cell", components, np.min(field), np.mean(field), np.max(field)]
        df = pd.DataFrame(data_array_dict).T
        return df.rename({0:"type", 1: "components", 2: "Min", 3: "Mean", 4: "Max"}, axis='columns')

    @property
    def cell_center_points(self):
        """
        Method for obtaining cell center points
        """
        if self._cell_center_points is not None:
            return self._cell_center_points
        ccf = vtk.vtkCellCenters()
        ccf.SetInputData(self.output)
        ccf.VertexCellsOn()
        ccf.Update()
        self._cell_center_points = vtk_to_numpy(ccf.GetOutput().GetPoints().GetData())
        if self.dim == 1:
            self.one_d_axis = self.one_d_axis
            self._cell_center_points = self._cell_center_points[:, self.one_d_axis]
        if self.dim == 2:
            self.plane = [0, 1, 2]
            self.plane.pop(self.two_d_planenormal)
            self._cell_center_points = np.delete(self._cell_center_points, self.two_d_planenormal, 1)
        return self._cell_center_points

    def get_neighbors(self, points_interpol, data_type="point"):
        """
        Method for obtaining neighbor points for interpolation.
        """
        points = self.points if data_type == "point" else self.cell_center_points
        df = pd.DataFrame(points)
        neighbors = {}
        if self.dim == 1:
            return neighbors
        for i, (_, val) in enumerate(points_interpol.items()):
            if self.dim == 2:
                x, y = self.plane
                df["r_"+str(i)] = (df[x]-val[x]) * (df[x]-val[x]) + (df[y]-val[y]) * (df[y]-val[y])
            elif self.dim == 3:
                df["r_"+str(i)] = ((df[0]-val[0]) * (df[0]-val[0]) + (df[1]-val[1]) * (df[1]-val[1])
                        + (df[2]-val[2]) * (df[2]-val[2]))
            neighbors[i] = df.sort_values(by=["r_" + str(i)]).head(self.nneighbors).index
        return neighbors

    def get_nearest_points(self, points_interpol):
        """
        Return a dictionary with closest mesh points

        Parameters
        ----------
        points_interpol : `dict`
        """
        if isinstance(points_interpol, dict):
            nb = self.get_neighbors(points_interpol)
            nearest = {}
            for i, (key, _) in enumerate(points_interpol.items()):
                nearest[key] = self.points[nb[i][0]]
        elif isinstance(points_interpol, np.ndarray):
            locator = vtk.vtkStaticPointLocator()
            locator.SetDataSet(self.output)
            locator.BuildLocator()
            nearestindices = []
            for point in points_interpol:
                nearestindices.append(locator.FindClosestPoint([point[0], point[1], point[2]]))
            nearest[nearestindices]
        return nearest

    def get_nearest_indices(self, points_interpol):
        """
        Return a dictionary with closest mesh point indices

        Parameters
        ----------
        points_interpol : `dict`
        """
        if isinstance(points_interpol, dict):
            nb = self.get_neighbors(points_interpol)
            nearest = {}
            for i, (key, _) in enumerate(points_interpol.items()):
                nearest[key] = nb[i][0]
        elif isinstance(points_interpol, np.ndarray):
            locator = vtk.vtkStaticPointLocator()
            locator.SetDataSet(self.output)
            locator.BuildLocator()
            nearest = []
            for point in points_interpol:
                nearest.append(locator.FindClosestPoint([point[0], point[1], point[2]]))
        return nearest


    def get_data_scipy(self, neighbors, points_interpol, fieldname, data_type="point", interpolation_method="linear"):
        """
        Get interpolated data for points_interpol using neighbor points.
        """
        field = self.get_point_field(fieldname) if data_type == "point" else self.get_cell_field(fieldname)
        points = self.points if data_type == "point" else self.cell_center_points
        resp = {}
        for i, (key, val) in enumerate(points_interpol.items()):
            if self.dim == 1:
                data = pd.DataFrame(points, columns = ['x'])
                data["y"] = field
                data.sort_values(by = ['x'], inplace=True)
                data.drop_duplicates(subset=['x'], inplace=True)
                f = interp1d(data['x'], data['y'])
                resp[key] = f(val[self.one_d_axis])
            elif self.dim == 2:
                x, y = self.plane
                grid_x, grid_y = np.array([[[val[x]]],[[val[y]]]])
                resp[key] = griddata(points[neighbors[i]], field[neighbors[i]],
                        (grid_x, grid_y), method=interpolation_method)[0][0]
            else:
                grid_x, grid_y, grid_z = np.array([[[[val[0]]]], [[[val[1]]]], [[[val[2]]]]])
                resp[key] = griddata(points[neighbors[i]], field[neighbors[i]],
                        (grid_x, grid_y, grid_z), method=interpolation_method)[0][0][0]
        return resp

    def get_data_vtk(self, points_interpol, interpolation_method="linear"):
        """
        Get interpolated data for points_interpol using vtks built-in interpolation methods
        """
        kernels = {"voronoi": vtk.vtkVoronoiKernel(), "gaussian": vtk.vtkGaussianKernel(),
        "shepard": vtk.vtkShepardKernel(), "linear": vtk.vtkLinearKernel()}
        pointnumpyarray = np.array([points_interpol[pt] for pt in points_interpol])
        out_u_grid = vtk.vtkUnstructuredGrid()
        r = vtk.vtkPoints()
        r.SetData(numpy_to_vtk(pointnumpyarray))
        out_u_grid.SetPoints(r)
        locator = vtk.vtkStaticPointLocator()
        locator.SetDataSet(self.output)
        locator.BuildLocator()
        interpolator = vtk.vtkPointInterpolator()
        interpolator.SetInputData(out_u_grid)
        interpolator.SetSourceData(self.output)
        interpolator.SetKernel(kernels[interpolation_method])
        interpolator.SetLocator(locator)
        interpolator.Update()
        return interpolator.GetOutput().GetPointData()



    def get_point_field(self, fieldname):
        """
        Return vtu cell field as numpy array.

        Parameters
        ----------
        fieldname : `str`
        """
        field = vtk_to_numpy(self.pdata.GetArray(fieldname))
        return field

    def get_cell_field(self, fieldname):
        """
        Return vtu point field as numpy array.

        Parameters
        ----------
        fieldname : `str`
        """
        field = vtk_to_numpy(self.cdata.GetArray(fieldname))
        return field

    def get_cell_field_as_point_data(self, fieldname):
        """
        Return vtu cell field as point field.

        Parameters
        ----------
        fieldname : `str`
        """
        c2p = vtk.vtkCellDataToPointData()
        c2p.SetInputData(self.output)
        c2p.Update()
        outpoints = c2p.GetOutput()
        nodes = outpoints.GetPointData()
        array =  vtk_to_numpy(nodes.GetArray(fieldname))
        return array

    def get_cell_field_names(self):
        """
        Get names of all cell fields in the vtu file.
        """
        fieldnames = []
        for i in range(self.cdata.GetNumberOfArrays()):
            fieldnames.append(self.cdata.GetArrayName(i))
        return fieldnames

    def get_point_field_names(self):
        """
        Get names of all point fields in the vtu file.
        """
        fieldnames = []
        for i in range(self.pdata.GetNumberOfArrays()):
            fieldnames.append(self.pdata.GetArrayName(i))
        return fieldnames

    def get_data(self, fieldname, pts = None, data_type="point", interpolation_method="linear"):
        """
        Get data of field "fieldname" at all points specified in "pts".

        Parameters
        ----------
        fieldname : `str`
        pts : `dict`, optional
              default: {'pt0': (0.0,0.0,0.0)}
        interpolation_method : `str`, optional
                               default: 'linear'
        """
        if pts is None:
            pts = {'pt0': (0.0,0.0,0.0)}
        resp = {}
        for pt in pts:
            if isinstance(fieldname, str):
                resp[pt] = []
            elif isinstance(fieldname, list):
                resp[pt] = {}
                for field in fieldname:
                    resp[pt][field] = []
        # TODO: move following part into separate method (similar code in PVDIO)
        if self.interpolation_backend == "scipy":
            nb = self.get_neighbors(pts, data_type=data_type)
            if isinstance(fieldname, str):
                data = self.get_data_scipy(nb, pts, fieldname, data_type=data_type,
                        interpolation_method=interpolation_method)
                for pt in pts:
                    resp[pt] = data[pt]
            elif isinstance(fieldname, list):
                data = {}
                for field in fieldname:
                    data[field] = self.get_data_scipy(nb, pts, field, data_type=data_type,
                            interpolation_method=interpolation_method)
                for pt in pts:
                    for field in fieldname:
                        resp[pt][field] = data[field][pt]
        elif self.interpolation_backend == "vtk":
            if data_type != "point":
                raise RuntimeError("reading cell data is not working with vtk backend yet")
            if isinstance(fieldname, str):
                resp_array = vtk_to_numpy(self.get_data_vtk(
                        pts, interpolation_method=interpolation_method).GetArray(fieldname))
                for i, pt in enumerate(pts):
                    resp[pt] = resp_array[i]

            elif isinstance(fieldname, list):
                resp_array_dict = {}
                vtkdata = self.get_data_vtk(pts, interpolation_method=interpolation_method)
                for field in fieldname:
                    resp_array_dict[field] = vtk_to_numpy(vtkdata.GetArray(fieldname))
                for i, pt in enumerate(pts):
                    for field in fieldname:
                        resp[pt][field] = resp_array_dict[field][i]
        else:
            raise RuntimeError(f"Interpolation backend {self.interpolation_backend} not found.")
        return resp


    def get_set_data(self, fieldname, pointsetarray=None, data_type="point", interpolation_method="linear"):
        """
        Get data specified in fieldname at all points specified in "pointsetarray".

        Parameters
        ----------
        fieldname : `str`
        pointsetarray : `list`, `numpy.ndarray` or `str`
                        if `str`, pointsetarray is construed as a filename of a submesh
        interpolation_method : `str`, optional
                               default: 'linear'
        """
        if pointsetarray is None:
            raise RuntimeError("No pointsetarray given.")
        if isinstance(pointsetarray, str):
            vtu = VTUIO(pointsetarray, dim=3)
            pointsetarray = vtu.points
        pts = {}
        # convert into point dictionary
        for i, entry in enumerate(pointsetarray):
            pts['pt'+str(i)] = entry
        resp = self.get_data(fieldname, pts=pts, data_type=data_type, interpolation_method=interpolation_method)
        resp_list = []
        # convert point dictionary into list
        for i, entry in enumerate(pointsetarray):
            resp_list.append(resp['pt' + str(i)])
        resp_array = np.array(resp_list)
        return resp_array

    def func_to_field(self, function, fieldname, ofilename, cell=False, writefile=True):
        """
        Add a field to the vtu file (which will be saved directly as "ofilename"
        by providing a three argument function(x,y,z)

        Parameters
        ----------
        function : `function`
        fieldname : `str`
        ofilename : `str`
        cell : `bool`
        """
        if callable(function) is False:
            print("function is not a function")
            raise TypeError
        if cell is True:
            points = self.cell_center_points
        else:
            points = self.points
        fieldarray = np.zeros(len(points))
        for i,_ in enumerate(fieldarray):
            if self.dim == 2:
                fieldarray[i] = function(points[i,0], points[i,1], 0.0)
            else:
                fieldarray[i] = function(points[i,0], points[i,1], points[i,2])
        field_vtk = numpy_to_vtk(fieldarray)
        if cell is True:
            r = self.cdata.AddArray(field_vtk)
            self.cdata.GetArray(r).SetName(fieldname)
        else:
            r = self.pdata.AddArray(field_vtk)
            self.pdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(ofilename)
            writer.SetInputData(self.output)
            writer.Write()

    def func_to_m_dim_field(self, functionarray, fieldname, ofilename, cell=False, writefile=True):
        """
        Add a multidimensional field to the vtu file (which will be saved directly as "ofilename"
        by providing am array of three argument functions.

        Parameters
        ----------
        functionarray : `array` of objects
        fieldname : `str`
        ofilename : `str`
        """
        mdim = len(functionarray)
        for function in functionarray:
            if callable(function) is False:
                print("functionarray is not containing functions only.")
                raise TypeError
        if cell is True:
            points = self.cell_center_points
        else:
            points = self.points
        fieldarray = np.zeros((len(points), mdim))
        for i,_ in enumerate(fieldarray):
            for j, func in enumerate(functionarray):
                if self.dim == 2:
                    fieldarray[i,j] = func(
                        points[i,0],
                        points[i,1],
                        0.0)
                else:
                    fieldarray[i,j] = func(
                        points[i,0],
                        points[i,1],
                        points[i,2])
        field_vtk = numpy_to_vtk(fieldarray)
        if cell is True:
            r = self.cdata.AddArray(field_vtk)
            self.cdata.GetArray(r).SetName(fieldname)
        else:
            r = self.pdata.AddArray(field_vtk)
            self.pdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(ofilename)
            writer.SetInputData(self.output)
            writer.Write()

    def point_data_to_cell_data(self, fieldname, ofilename, writefile=True):
        """
        convert pointdata to cell data of field "fieldname"

        Parameters
        ----------
        fieldname : `str`
        ofilename : `str`
        """
        p2c = vtk.vtkPointDataToCellData()
        p2c.SetInputData(self.output)
        p2c.Update()
        outcells = p2c.GetOutput()
        cells = outcells.GetCellData()
        array =  cells.GetArray(fieldname)
        cells_orig = self.output.GetCellData()
        cells_orig.AddArray(array)
        if writefile is True:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(ofilename)
            writer.SetInputData(self.output)
            writer.Write()

    def delete_point_field(self, fieldnames, ofilename, writefile=True):
        """
        delete point field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
        ofilename : `str`
        """
        if isinstance(fieldnames, str):
            self.pdata.RemoveArray(fieldnames)
        elif isinstance(fieldnames, list):
            for fieldname in fieldnames:
                self.pdata.RemoveArray(fieldname)
        else:
            raise TypeError("Fieldnames has the wrong type. Please provide a list or string.")
        if writefile is True:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(ofilename)
            writer.SetInputData(self.output)
            writer.Write()

    def delete_cell_field(self, fieldnames, ofilename, writefile=True):
        """
        delete cell field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
        ofilename : `str`
        """
        if isinstance(fieldnames, str):
            self.cdata.RemoveArray(fieldnames)
        elif isinstance(fieldnames, list):
            for fieldname in fieldnames:
                self.cdata.RemoveArray(fieldname)
        else:
            raise TypeError("Fieldnames has the wrong type. Please provide a list or string.")
        if writefile is True:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(ofilename)
            writer.SetInputData(self.output)
            writer.Write()

    def write_point_field(self, field, fieldname, ofilename, writefile=True):
        """
        Write a field (numpy array of correct size)
        to field "fieldname" as file "ofilename".

        Parameters
        ----------
        field : `array`
        fieldname : `str`
        ofilename : `str`
        """
        field_vtk = numpy_to_vtk(field)
        r = self.pdata.AddArray(field_vtk)
        self.pdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(ofilename)
            writer.SetInputData(self.output)
            writer.Write()

    def write_cell_field(self, field, fieldname, ofilename, writefile=True):
        """
        Write a field (numpy array of correct size)
        to field "fieldname" as file "ofilename".

        Parameters
        ----------
        field : `array`
        fieldname : `str`
        ofilename : `str`
        """
        field_vtk = numpy_to_vtk(field)
        r = self.cdata.AddArray(field_vtk)
        self.cdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(ofilename)
            writer.SetInputData(self.output)
            writer.Write()

    def write(self, filename):
        """
        Write data as file "filename".

        Parameters
        ----------
        filename : `str`
        """
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.output)
        writer.Write()



class PVDIO:
    """
    Class for handling I/O of PVD files

    Parameters
    ----------
    filename : `str`
    nneighbors : `int`, optional
    dim : `int`
    one_d_axis : `int`
                 between 0 and 2, default: 0
    two_d_planenormal : `int`
                  between 0 and 2, default: 2
    interpolation_backend : `str`
                 scipy or vtk
    """
    def __init__(self, filename, nneighbors=20, dim=3, one_d_axis=0, two_d_planenormal=2,
                                                            interpolation_backend="scipy"):
        if os.path.isfile(filename) is True:
            self.folder, self.filename = os.path.split(filename)
        else:
            raise RuntimeError(f"File not found: {filename}")
        self.nneighbors = nneighbors
        self.timesteps = np.array([])
        self.vtufilenames = []
        self.read_pvd(os.path.join(self.folder, self.filename))
        self.dim = dim
        self.one_d_axis = one_d_axis
        self.two_d_planenormal = two_d_planenormal
        self.interpolation_backend = interpolation_backend

    def delete_point_field(self, fieldnames):
        """
        delete point field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
        """
        for filename in self.vtufilenames:
            vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
            vtu.delete_point_field(fieldnames, filename)

    def delete_cell_field(self, fieldnames):
        """
        delete cell field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
        """
        for filename in self.vtufilenames:
            vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
            vtu.delete_cell_field(fieldnames, filename)

    def read_pvd(self, filename):
        """
        Read in PVD file

        Parameters
        ----------
        filename : `str`
        """
        self.filename = filename
        tree = ET.parse(self.filename)
        root = tree.getroot()
        for collection in root.getchildren():
            for dataset in collection.getchildren():
                self.timesteps = np.append(self.timesteps, [float(dataset.attrib['timestep'])])
                self.vtufilenames.append(dataset.attrib['file'])

    def read_time_series(self, fieldname, pts=None, data_type="point", interpolation_method="linear"):
        """
        Return time series data of field "fieldname" at points pts.
        Also a list of fieldnames can be provided as "fieldname"

        Parameters
        ----------
        fieldname : `str`
        pts : `dict`, optional
        data_type : `str` optional
              "point" or "cell"
        interpolation_method : `str`, optional
                               default: 'linear
        """
        if pts is None:
            raise RuntimeError("No points given")
        resp_t = {}
        for pt in pts:
            if isinstance(fieldname, str):
                resp_t[pt] = []
            elif isinstance(fieldname, list):
                resp_t[pt] = {}
                for field in fieldname:
                    resp_t[pt][field] = []
        for i, filename in enumerate(self.vtufilenames):
            vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
            if self.interpolation_backend == "scipy":
                if i == 0:
                    nb = vtu.get_neighbors(pts, data_type=data_type)
                if isinstance(fieldname, str):
                    data = vtu.get_data_scipy(nb, pts, fieldname, data_type=data_type,
                            interpolation_method=interpolation_method)
                    for pt in pts:
                        resp_t[pt].append(data[pt])
                elif isinstance(fieldname, list):
                    data = {}
                    for field in fieldname:
                        data[field] = vtu.get_data_scipy(nb, pts, field, data_type=data_type,
                                interpolation_method=interpolation_method)
                    for pt in pts:
                        for field in fieldname:
                            resp_t[pt][field].append(data[field][pt])
            elif self.interpolation_backend == "vtk":
                if data_type != "point":
                    raise RuntimeError("reading cell data is not working with vtk backend yet")
                if isinstance(fieldname, str):
                    data = vtk_to_numpy(
                        vtu.get_data_vtk(pts, interpolation_method=interpolation_method).GetArray(fieldname))
                    for j, pt in enumerate(pts):
                        resp_t[pt].append(data[j])
                elif isinstance(fieldname, list):
                    data = {}
                    vtkdata = vtu.get_data_vtk(pts, interpolation_method=interpolation_method)
                    for field in fieldname:
                        data[field] = vtk_to_numpy(vtkdata.GetArray(fieldname))
                    for j, pt in enumerate(pts):
                        for field in fieldname:
                            resp_t[pt][field].append(data[field][j])
        resp_t_array = {}
        for pt, field in resp_t.items():
            if isinstance(fieldname, str):
                resp_t_array[pt] = np.array(field)
            elif isinstance(fieldname, list):
                resp_t_array[pt] = {}
                for field_, fieldarray in resp_t[pt].items():
                    resp_t_array[pt][field_] = np.array(fieldarray)
        return resp_t_array

    def read_time_slice(self, time, fieldname):
        """
        Print field "fieldname" at time "time".

        Parameters
        ----------
        time : `float`
        fieldname : `str`
        """
        filename = None
        for i, ts in enumerate(self.timesteps):
            if time == ts:
                filename = self.vtufilenames[i]
        if not filename is None:
            vtu = VTUIO(os.path.join(self.folder,filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
            field = vtu.get_point_field(fieldname)
        else:
            filename1 = None
            filename2 = None
            time1 = 0.0
            time2 = 0.0
            for i, ts in enumerate(self.timesteps):
                try:
                    if ts < time < self.timesteps[i+1]:
                        time1 = ts
                        time2 = self.timesteps[i+1]
                        filename1 = self.vtufilenames[i]
                        filename2 = self.vtufilenames[i+1]
                except IndexError:
                    print("time is out of range")
            if (filename1 is None) or (filename2 is None):
                print("time is out of range")
            else:
                vtu1 = VTUIO(os.path.join(self.folder,filename1),
                        nneighbors=self.nneighbors, dim=self.dim,
                        one_d_axis=self.one_d_axis,
                        two_d_planenormal=self.two_d_planenormal,
                        interpolation_backend=self.interpolation_backend)
                vtu2 = VTUIO(os.path.join(self.folder,filename2),
                        nneighbors=self.nneighbors, dim=self.dim,
                        one_d_axis=self.one_d_axis,
                        two_d_planenormal=self.two_d_planenormal,
                        interpolation_backend=self.interpolation_backend)
                field1 = vtu1.get_point_field(fieldname)
                field2 = vtu2.get_point_field(fieldname)
                fieldslope = (field2-field1)/(time2-time1)
                field = field1 + fieldslope * (time-time1)
        return field

    def read_set_data(self, time, fieldname, pointsetarray = None, data_type="point", interpolation_method="linear"):
        """
        Get data of field "fieldname" at time "time" alon a given "pointsetarray".

        Parameters
        ----------
        time : `float`
        fieldname : `str`
        pointsetarray : `list`, `numpy.ndarray` or `str`
                        if `str`, pointsetarray is construed as a filename of a submesh
        interpolation_method : `str`
                               default: 'linear'
        """
        if pointsetarray is None:
            raise RuntimeError("No pointsetarray given.")
        if isinstance(pointsetarray, str):
            vtu = VTUIO(pointsetarray, dim=3)
            pointsetarray = vtu.points
        filename = None
        for i, ts in enumerate(self.timesteps):
            if time == ts:
                filename = self.vtufilenames[i]
        if not filename is None:
            vtu = VTUIO(os.path.join(self.folder,filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
            field = vtu.get_set_data(fieldname, pointsetarray, data_type=data_type, interpolation_method=interpolation_method)
        else:
            filename1 = None
            filename2 = None
            time1 = 0.0
            time2 = 0.0
            for i, ts in enumerate(self.timesteps):
                try:
                    if ts < time < self.timesteps[i+1]:
                        time1 = ts
                        time2 = self.timesteps[i+1]
                        filename1 = self.vtufilenames[i]
                        filename2 = self.vtufilenames[i+1]
                except IndexError:
                    print("time is out of range")
            if (filename1 is None) or (filename2 is None):
                print("time is out of range")
            else:
                vtu1 = VTUIO(os.path.join(self.folder,filename1),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
                vtu2 = VTUIO(os.path.join(self.folder,filename2),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
                field1 = vtu1.get_set_data(fieldname, pointsetarray, data_type=data_type, interpolation_method=interpolation_method)
                field2 = vtu2.get_set_data(fieldname, pointsetarray, data_type=data_type, interpolation_method=interpolation_method)
                fieldslope = (field2-field1)/(time2-time1)
                field = field1 + fieldslope * (time-time1)
        return field

    def read_aggregate(self, fieldname, agg_fct, data_type="point", pointsetarray=None):
        """
        Return time series data of an aggregate function for field "fieldname".

        Parameters
        ----------
        fieldname : `str` or `list`
        agg_fct : `str`,
              can be: "min", "max" or "mean"
        data_type : `str` optional
              "point" or "cell"
        pointsetarray : `str`, `list` or `numpy.ndarray`
                        defines a submesh
                        if `str` pointsetarray is construed as filename containing the mesh
        """
        agg_fcts = {"min": np.min,
                    "max": np.max,
                    "mean": np.mean}
        resp_t = {}
        if isinstance(fieldname, str):
            resp_t = []
        elif isinstance(fieldname, list):
            resp_t = {}
            for field in fieldname:
                resp_t[field] = []
        if not pointsetarray is None:
            if isinstance(pointsetarray, str):
                vtu = VTUIO(pointsetarray, dim=3)
                pointsetarray = vtu.points
            pointsetarray = np.array(pointsetarray)
        submeshindices = None
        for i, filename in enumerate(self.vtufilenames):
            vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
            if (i == 0) and (not pointsetarray is None):
                submeshindices = vtu.get_nearest_indices(pointsetarray)
            if isinstance(fieldname, str):
                if data_type == "point":
                    data = agg_fcts[agg_fct](vtu.get_point_field(fieldname)[submeshindices])
                elif data_type == "cell":
                    data = agg_fcts[agg_fct](vtu.get_cell_field(fieldname)[submeshindices])
                resp_t.append(data)
            elif isinstance(fieldname, list):
                for field in fieldname:
                    if data_type == "point":
                        data = agg_fcts[agg_fct](vtu.get_point_field(field)[submeshindices])
                    elif data_type == "cell":
                        data = agg_fcts[agg_fct](vtu.get_cell_field(field)[submeshindices])
                    resp_t[field].append(data)
        return resp_t

    def clear_pvd_rel_path(self, write=True):
        """
        Delete relative directory paths in the vtu filenames of the PVD file.

        Parameters
        ----------
        write : `bool`
        """
        xpath="./Collection/DataSet"
        tree = ET.parse(self.filename)
        root = tree.getroot()
        find_xpath = root.findall(xpath)
        for tag in find_xpath:
            filename = tag.get("file")
            filename_new = filename.split("/")[-1]
            tag.set("file", filename_new)
        if write is True:
            tree.write(self.filename,
                            encoding="ISO-8859-1",
                            xml_declaration=True,
                            pretty_print=True)
        #update file list:
        newlist = []
        for entry in  self.vtufilenames:
            newlist.append(entry.split("/")[-1])
        self.vtufilenames = newlist

    def write_xdmf(self, filename):
        import meshio
        mesh = meshio.read(self.vtufilenames[0])
        with meshio.xdmf.TimeSeriesWriter(filename) as writer:
            for i, t in enumerate(self.timesteps):
                mesh = meshio.read(self.vtufilenames[i])
                if i == 0:
                    writer.write_points_cells(mesh.points, mesh.cells)
                writer.write_data(t, point_data=mesh.point_data, cell_data=mesh.cell_data)
