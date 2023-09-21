# -*- coding: utf-8 -*-
"""
VTUinterface is a python package for easy accessing VTU/PVD files as
outputed by Finite Element software like OpenGeoSys. It uses the VTK python
wrapper and linear interpolation between time steps and grid points access
any points in and and time within the simulation domain.

Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
            Distributed under a Modified BSD License.
              See accompanying file LICENSE or
              http://www.opengeosys.org/project/license

"""

# pylint: disable=C0103, R0902, R0914, R0913
import os
import sys
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
                                                        interpolation_backend=None):
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
            self.ipdata = self.output.GetFieldData()
        except AttributeError:
            print(f"File {self.filename} does not contain any data")
        try:
            self.points = vtk_to_numpy(self.output.GetPoints().GetData())
        except AttributeError:
            print(f"File {self.filename} does not contain any points")
        self._cell_center_points = None
        self._plane = [0, 1, 2]
        self.dim = dim
        self.nneighbors = nneighbors
        self.one_d_axis=one_d_axis
        self.two_d_planenormal = two_d_planenormal
        self._plane = [0, 1, 2]
        self._plane.pop(self.two_d_planenormal)
        # interpolation settings
        if interpolation_backend is None:
            self.interpolation_backend = "vtk"
        else:
            self.interpolation_backend = interpolation_backend
        self.vtk_gaussian_sharpness = 4.0
        self.vtk_gaussian_radius = 0.5
        self.vtk_gaussian_footprint_to_n_closest = False
        self.vtk_shepard_power_parameter = 2.0
        self.vtk_shepard_radius = 0.5
        self.datamode = "appended"


    def _proj(self, array):
        if self.dim == 1:
            self.one_d_axis = self.one_d_axis
            array = array[:,self.one_d_axis]
        if self.dim == 2:
            array = np.delete(array, self.two_d_planenormal, 1)
        return array
    @property
    def header(self):
        header_dict = {"N Cells": [str(len(self.cell_center_points))], "N Points": [str(len(self.points))],
                "N Cell Arrays": [len(self.get_cell_field_names())],
                       "N Point Arrays": [len(self.get_point_field_names())], "N Integration Point Arrays": [len(self.get_integration_point_field_names())],
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
        ip_names = self.get_integration_point_field_names()
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
            if field.ndim == 1:
                components = 1
            else:
                components = field.shape[1]
            data_array_dict[name] = ["cell", components, np.min(field), np.mean(field), np.max(field)]
        for name in ip_names:
            field = self.get_integration_point_field(name)
            if field.ndim == 1:
                components = 1
            else:
                components = field.shape[1]
            data_array_dict[name] = ["ip data", components, np.min(field), np.mean(field), np.max(field)]
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
                x, y = self._plane
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
        points = self._proj(points)
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
                x, y = self._plane
                grid_x, grid_y = np.array([[[val[x]]],[[val[y]]]])
                resp[key] = griddata(points[neighbors[i]], field[neighbors[i]],
                        (grid_x, grid_y), method=interpolation_method)[0][0]
            else:
                grid_x, grid_y, grid_z = np.array([[[[val[0]]]], [[[val[1]]]], [[[val[2]]]]])
                resp[key] = griddata(points[neighbors[i]], field[neighbors[i]],
                        (grid_x, grid_y, grid_z), method=interpolation_method)[0][0][0]
        return resp

    def get_data_vtk(self, points_interpol, data_type="point", interpolation_method="probefilter"):
        """
        Get interpolated data for points_interpol using vtks built-in interpolation methods
        """
        kernels = {"nearest": vtk.vtkVoronoiKernel(), "voronoi": vtk.vtkVoronoiKernel(), "gaussian": vtk.vtkGaussianKernel(),
        "shepard": vtk.vtkShepardKernel()}
        #           "linear": vtk.vtkLinearKernel()} temporarily deactivated
        pointnumpyarray = np.array([points_interpol[pt] for pt in points_interpol])
        out_u_grid = vtk.vtkUnstructuredGrid()
        r = vtk.vtkPoints()
        r.SetData(numpy_to_vtk(pointnumpyarray))
        out_u_grid.SetPoints(r)
        locator = vtk.vtkStaticPointLocator()
        locator.SetDataSet(self.output)
        locator.BuildLocator()
        if interpolation_method == "probefilter":
            interpolator = vtk.vtkProbeFilter()
        else:
            interpolator = vtk.vtkPointInterpolator()
        interpolator.SetInputData(out_u_grid)
        if data_type == "point":
            interpolator.SetSourceData(self.output)
        else:
            out_u_grid_source = vtk.vtkUnstructuredGrid()
            r_source = vtk.vtkPoints()
            r_source.SetData(numpy_to_vtk(self.cell_center_points))
            out_u_grid_source.SetPoints(r_source)
            pdata_source = out_u_grid_source.GetPointData()
            cellfieldnames = self.get_cell_field_names()
            for fieldname in cellfieldnames:
                carray = self.cdata.GetArray(fieldname)
                q = pdata_source.AddArray(carray)
                pdata_source.GetArray(q).SetName(fieldname)
            interpolator.SetSourceData(out_u_grid_source)
        if interpolation_method != "probefilter":
            kernel = kernels[interpolation_method]
        if interpolation_method == "gaussian":
            if self.vtk_gaussian_footprint_to_n_closest is True:
                kernel.SetKernelFootprintToNClosest()
                kernel.SetNumberOfPoints(self.nneighbors)
            kernel.SetSharpness(self.vtk_gaussian_sharpness)
            kernel.SetRadius(self.vtk_gaussian_radius)
        if interpolation_method == "shepard":
            kernel.SetPowerParameter(self.vtk_shepard_power_parameter)
            kernel.SetRadius(self.vtk_shepard_radius)
        if not interpolation_method == "probefilter":
            interpolator.SetKernel(kernel)
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

    def get_integration_point_field(self, fieldname):
        """
        Return integration point field.

        Parameters
        ----------
        fieldname : `str`
        """
        field = vtk_to_numpy(self.ipdata.GetArray(fieldname))
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

    def get_integration_point_field_names(self):
        """
        Get names of all integration point data fields in the vtu file.
        """
        fieldnames = []
        for i in range(self.ipdata.GetNumberOfArrays()):
            name = self.ipdata.GetArrayName(i)
            if "_ip" in name:
                fieldnames.append(name)
        return fieldnames

    def get_data(self, fieldname, pts = None, data_type="point", interpolation_method=None):
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
            if interpolation_method is None:
                interpolation_method="linear"
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
            if interpolation_method is None:
                interpolation_method="probefilter"
            if isinstance(fieldname, str):
                resp_array = vtk_to_numpy(self.get_data_vtk(
                        pts, data_type=data_type, interpolation_method=interpolation_method).GetArray(fieldname))
                for i, pt in enumerate(pts):
                    resp[pt] = resp_array[i]

            elif isinstance(fieldname, list):
                resp_array_dict = {}
                vtkdata = self.get_data_vtk(pts, data_type=data_type, interpolation_method=interpolation_method)
                for field in fieldname:
                    resp_array_dict[field] = vtk_to_numpy(vtkdata.GetArray(field))
                for i, pt in enumerate(pts):
                    for field in fieldname:
                        resp[pt][field] = resp_array_dict[field][i]
        else:
            raise RuntimeError(f"Interpolation backend {self.interpolation_backend} not found.")
        return resp


    def get_set_data(self, fieldname, pointsetarray=None, data_type="point", interpolation_method=None):
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

    def func_to_field(self, function, fieldname, ofilename=None, cell=False, array_type=None, writefile=True, datamode=None):
        """
        Add a field to the vtu file (which will be saved directly as "ofilename"
        by providing a three argument function(x,y,z)

        Parameters
        ----------
        function : `function`
        fieldname : `str`
        ofilename : `str`
        cell : `bool`
        array_type : `vtk array type`
        writefile : `bool`
        datamode : `str`
        """
        if ofilename is None:
            ofilename = self.filename
        if callable(function) is False:
            print("function is not a function")
            raise TypeError
        if cell is True:
            points = self.cell_center_points
        else:
            points = self.points
        fieldarray = np.zeros(len(points))
        for i,_ in enumerate(fieldarray):
            if self.dim == 1:
                fieldarray[i] = function(points[i,0], 0.0, 0.0)
            elif self.dim == 2:
                fieldarray[i] = function(points[i,0], points[i,1], 0.0)
            else:
                fieldarray[i] = function(points[i,0], points[i,1], points[i,2])
        field_vtk = numpy_to_vtk(fieldarray, array_type=array_type)
        if cell is True:
            r = self.cdata.AddArray(field_vtk)
            self.cdata.GetArray(r).SetName(fieldname)
        else:
            r = self.pdata.AddArray(field_vtk)
            self.pdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            self.write(ofilename, datamode=datamode)

    def func_to_m_dim_field(self, functionarray, fieldname, ofilename=None, cell=False, array_type=None, writefile=True, datamode=None):
        """
        Add a multidimensional field to the vtu file (which will be saved directly as "ofilename"
        by providing am array of three argument functions.

        Parameters
        ----------
        functionarray : `array` of objects
        fieldname : `str`
        ofilename : `str`
        cell : `bool`
        array_type : `vtk array type`
        writefile : `bool`
        data_mode : `str`
        """
        if ofilename is None:
            ofilename = self.filename
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
                if self.dim == 1:
                    fieldarray[i,j] = func(
                        points[i,0],
                        0.0,
                        0.0)
                elif self.dim == 2:
                    fieldarray[i,j] = func(
                        points[i,0],
                        points[i,1],
                        0.0)
                else:
                    fieldarray[i,j] = func(
                        points[i,0],
                        points[i,1],
                        points[i,2])
        field_vtk = numpy_to_vtk(fieldarray, array_type=array_type)
        if cell is True:
            r = self.cdata.AddArray(field_vtk)
            self.cdata.GetArray(r).SetName(fieldname)
        else:
            r = self.pdata.AddArray(field_vtk)
            self.pdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            self.write(ofilename, datamode=datamode)

    def point_data_to_cell_data(self, fieldname, ofilename=None, writefile=True, datamode=None):
        """
        convert pointdata to cell data of field "fieldname"

        Parameters
        ----------
        fieldname : `str`
        ofilename : `str`
        datamode : `str`
        """
        if ofilename is None:
            ofilename = self.filename
        p2c = vtk.vtkPointDataToCellData()
        p2c.SetInputData(self.output)
        p2c.Update()
        outcells = p2c.GetOutput()
        cells = outcells.GetCellData()
        array =  cells.GetArray(fieldname)
        cells_orig = self.output.GetCellData()
        cells_orig.AddArray(array)
        if writefile is True:
            self.write(ofilename, datamode=datamode)

    def delete_point_field(self, fieldnames=None, ofilename=None, writefile=True, datamode=None):
        """
        delete point field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
            if `None` all fields will be deleted
        ofilename : `str`
        writefile : `bool`
        datamode : `str`
        """
        if fieldnames is None:
            fieldnames = self.get_point_field_names()
        if ofilename is None:
            ofilename = self.filename
        if isinstance(fieldnames, str):
            self.pdata.RemoveArray(fieldnames)
        elif isinstance(fieldnames, list):
            for fieldname in fieldnames:
                self.pdata.RemoveArray(fieldname)
        else:
            raise TypeError("Fieldnames has the wrong type. Please provide a list or string.")
        if writefile is True:
            self.write(ofilename, datamode=datamode)

    def delete_cell_field(self, fieldnames=None, ofilename=None, writefile=True, datamode=None):
        """
        delete cell field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
                if `None` all fields will be deleted
        ofilename : `str`
        writefile : `bool`
        datamode : `str`
        """
        if fieldnames is None:
            fieldnames = self.get_cell_field_names()
        if ofilename is None:
            ofilename = self.filename
        if isinstance(fieldnames, str):
            self.cdata.RemoveArray(fieldnames)
        elif isinstance(fieldnames, list):
            for fieldname in fieldnames:
                self.cdata.RemoveArray(fieldname)
        else:
            raise TypeError("Fieldnames has the wrong type. Please provide a list or string.")
        if writefile is True:
            self.write(ofilename, datamode=datamode)

    def delete_integration_point_field(self, fieldnames=None, ofilename=None, writefile=True, datamode=None):
        """
        delete integration point field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
            if `None` all fields will be deleted
        ofilename : `str`
        writefile : `bool`
        datamode : `str`
        """
        if fieldnames is None:
            fieldnames = self.get_integration_point_field_names()
        if ofilename is None:
            ofilename = self.filename
        if isinstance(fieldnames, str):
            self.ipdata.RemoveArray(fieldnames)
        elif isinstance(fieldnames, list):
            for fieldname in fieldnames:
                self.ipdata.RemoveArray(fieldname)
        else:
            raise TypeError("Fieldnames has the wrong type. Please provide a list or string.")
        if writefile is True:
            self.write(ofilename, datamode=datamode)


    def add_point_field(self, field, fieldname, ofilename, writefile=True, array_type=None, datamode=None):
        """
        Write a field (numpy array of correct size)
        to field "fieldname" as file "ofilename".

        Parameters
        ----------
        field : `array`
        fieldname : `str`
        ofilename : `str`
        writefile : `bool`
        array_type : `vtk array type`
        datamode : `str`
        """
        field_vtk = numpy_to_vtk(field, array_type=array_type)
        r = self.pdata.AddArray(field_vtk)
        self.pdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            self.write(ofilename, datamode=datamode)

    def add_cell_field(self, field, fieldname, ofilename, writefile=True, array_type=None, datamode=None):
        """
        Write a field (numpy array of correct size)
        to field "fieldname" as file "ofilename".

        Parameters
        ----------
        field : `array`
        fieldname : `str`
        ofilename : `str`
        writefile : `bool`
        array_type : `vtk array type`
        datamode : `str`
        """
        field_vtk = numpy_to_vtk(field)
        r = self.cdata.AddArray(field_vtk)
        self.cdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            self.write(ofilename, datamode=datamode)

    def add_integration_point_field(self, field, fieldname, ofilename, writefile=True, array_type=None, datamode=None):
        """
        Write a field (numpy array of correct size)
        to field "fieldname" as file "ofilename".

        Parameters
        ----------
        field : `array`
        fieldname : `str`
        ofilename : `str`
        writefile : `bool`
        array_type : `vtk array type`
        datamode : `str`
        """
        field_vtk = numpy_to_vtk(field)
        r = self.ipdata.AddArray(field_vtk)
        self.ipdata.GetArray(r).SetName(fieldname)
        if writefile is True:
            self.write(ofilename, datamode=datamode)

    def write(self, filename, datamode=None):
        """
        Write data as file "filename".

        Parameters
        ----------
        filename : `str`
        datamode : `str`
        """
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.output)
        if not datamode is None:
            self.datamode = datamode
        if self.datamode == "binary":
            writer.SetDataModeToBinary()
        elif self.datamode == "ascii":
            writer.SetDataModeToAscii()
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
                                                            interpolation_backend=None):
        if interpolation_backend is None:
            self.interpolation_backend = "vtk"
        else:
            self.interpolation_backend = interpolation_backend
        if os.path.isfile(filename) is True:
            self.folder, self.filename = os.path.split(filename)
        else:
            raise RuntimeError(f"File not found: {filename}")
        self.nneighbors = nneighbors
        self.timesteps = np.array([])
        self.vtufilenames = []
        self.tree = None
        self.read_pvd(os.path.join(self.folder, self.filename))
        self.dim = dim
        self.one_d_axis = one_d_axis
        self.two_d_planenormal = two_d_planenormal

    def get_cell_field_names(self):
        """
        Get names of all cell fields in the vtu file.
        """
        filename = self.vtufilenames[0]
        vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
        return vtu.get_cell_field_names()

    def get_point_field_names(self):
        """
        Get names of all point fields in the vtu file.
        """
        filename = self.vtufilenames[0]
        vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
        return vtu.get_point_field_names()

    def get_integration_point_field_names(self):
        """
        Get names of integration point fields in the vtu file.
        """
        filename = self.vtufilenames[0]
        vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
        return vtu.get_integration_point_field_names()

    def delete_point_field(self, fieldnames=None):
        """
        delete point field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
            if `None` all fields will be deleted
        """
        for filename in self.vtufilenames:
            vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
            vtu.delete_point_field(fieldnames, filename)

    def delete_cell_field(self, fieldnames=None):
        """
        delete cell field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
            if `None` all fields will be deleted
        """
        for filename in self.vtufilenames:
            vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
            vtu.delete_cell_field(fieldnames, filename)

    def delete_integration_point_field(self, fieldnames=None, skip_last=False):
        """
        delete integration point field(s) and write data to disk

        Parameters
        ----------
        fieldnames : `str` or `list`
                if `None` all fields will be deleted
        skip_last : `boolean`
        """
        nmax = len(self.vtufilenames)
        for i, filename in enumerate(self.vtufilenames):
            if ((i+skip_last) < nmax):
                vtu = VTUIO(os.path.join(self.folder, filename),
                    nneighbors=self.nneighbors, dim=self.dim,
                    one_d_axis=self.one_d_axis,
                    two_d_planenormal=self.two_d_planenormal,
                    interpolation_backend=self.interpolation_backend)
                vtu.delete_integration_point_field(fieldnames, os.path.join(self.folder,filename))

    def read_pvd(self, filename):
        """
        Read in PVD file

        Parameters
        ----------
        filename : `str`
        """
        self.tree = ET.parse(filename)
        root = self.tree.getroot()
        for collection in root.getchildren():
            for dataset in collection.getchildren():
                self.timesteps = np.append(self.timesteps, [float(dataset.attrib['timestep'])])
                self.vtufilenames.append(dataset.attrib['file'])

    def read_time_series(self, fieldname, pts=None, data_type="point", interpolation_method=None):
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
                if interpolation_method is None:
                    interpolation_method = "linear"
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
                if interpolation_method is None:
                      interpolation_method = "probefilter"
                if isinstance(fieldname, str):
                    data = vtk_to_numpy(
                        vtu.get_data_vtk(pts, data_type=data_type, interpolation_method=interpolation_method).GetArray(fieldname))
                    for j, pt in enumerate(pts):
                        resp_t[pt].append(data[j])
                elif isinstance(fieldname, list):
                    data = {}
                    vtkdata = vtu.get_data_vtk(pts, data_type=data_type, interpolation_method=interpolation_method)
                    for field in fieldname:
                        data[field] = vtk_to_numpy(vtkdata.GetArray(field))
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

    def read_set_data(self, time, fieldname, pointsetarray = None, data_type="point", interpolation_method=None):
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
            return field
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

    def clear_pvd_path(self, write=True):
        """
        Delete relative and absolute directory paths in the vtu filenames of the PVD file.

        Parameters
        ----------
        write : `bool`
        """
        xpath="./Collection/DataSet"
        tree = ET.parse(os.path.join(self.folder, self.filename))
        root = tree.getroot()
        find_xpath = root.findall(xpath)
        for tag in find_xpath:
            filename = tag.get("file")
            filename_new = filename.split("/")[-1]
            tag.set("file", filename_new)
        if write is True:
            tree.write(os.path.join(self.folder, self.filename),
                            encoding="ISO-8859-1",
                            xml_declaration=True,
                            pretty_print=True)
        #update file list:
        newlist = []
        for entry in  self.vtufilenames:
            if sys.platform == "win32":
                newlist.append(entry.split("\\")[-1])
            else:
                # TODO: Check function on other operating systems
                newlist.append(entry.split("/")[-1])
        self.vtufilenames = newlist

    def rename(self, newname):
        """
        Rename PVD file

        Parameters
        ----------
        newname : `str`
        """
        tree = ET.parse(os.path.join(self.folder, self.filename))
        if not ".pvd" in newname:
            newname = newname + ".pvd"
        os.rename(os.path.join(self.folder, self.filename), os.path.join(self.folder, newname))
        newname = newname.split(".pvd")[0]
        vtufilelist = self.vtufilenames
        for i, filename in enumerate(vtufilelist):
            newvtuname = filename.replace(self.filename.split(".pvd")[0], newname)
            self.vtufilenames[i] = newvtuname
            os.rename(os.path.join(self.folder, filename), os.path.join(self.folder, newvtuname))
        xpath="./Collection/DataSet"
        root = tree.getroot()
        find_xpath = root.findall(xpath)
        for tag in find_xpath:
            filename = tag.get("file")
            filename_new = filename.replace(self.filename.split(".pvd")[0], newname)
            tag.set("file", filename_new)
        self.filename = f"{newname}.pvd"
        tree.write(os.path.join(self.folder, self.filename),
                            encoding="ISO-8859-1",
                            xml_declaration=True,
                            pretty_print=True)

    def append(self, filename, vtu_rename=False):
        """
        appends entries from another PVD file

        Parameters
        ----------
        filename : `str`
        vtu_rename : `bool`
        """
        tree = ET.parse(filename)
        xpath="./Collection/DataSet"
        root = tree.getroot()
        find_xpath = root.findall(xpath)
        offset = 0
        folder, filename_dirstripped = os.path.split(filename)
        if float(find_xpath[0].get("timestep")) < self.timesteps[-1]:
            offset = self.timesteps[-1]
        elif float(find_xpath[0].get("timestep")) == self.timesteps[-1]:
            self.timesteps = self.timesteps[:-1]
            self.vtufilenames = self.vtufilenames[:-1]
        for tag in find_xpath:
            self.timesteps = np.append(self.timesteps, [float(tag.attrib['timestep'])+offset])
            newvtuname = tag.attrib['file']
            if vtu_rename is True:
                newvtuname = tag.attrib['file'].replace(filename_dirstripped.split(".pvd")[0],
                                                   self.filename.split(".pvd")[0])
                os.rename(os.path.join(self.folder, tag.attrib['file']), os.path.join(self.folder,
                                                                     folder, newvtuname))
            self.vtufilenames.append(os.path.join(self.folder, folder, newvtuname))
        root = ET.Element("VTKFile")
        root.attrib["type"] = "Collection"
        root.attrib["version"] = "0.1"
        root.attrib["byte_order"] = "LittleEndian"
        root.attrib["compressor"] = "vtkZLibDataCompressor"
        collection = ET.SubElement(root,"Collection")
        timestepselements = []
        #pvdwriter
        for i, timestep in enumerate(self.timesteps):
            timestepselements.append(ET.SubElement(collection, "DataSet"))
            timestepselements[-1].attrib["timestep"] = str(timestep)
            timestepselements[-1].attrib["group"] = ""
            timestepselements[-1].attrib["part"] = "0"
            timestepselements[-1].attrib["file"] = self.vtufilenames[i]
        tree = ET.ElementTree(root)
        tree.write(os.path.join(self.folder, self.filename), encoding="ISO-8859-1",
                   xml_declaration=True, pretty_print=True)

    @staticmethod
    def create_pvd_from_pattern(filename, pattern=""):
        files = os.listdir()
        vtufiles = []
        timesteps = []
        for file in files:
            if (pattern in file) and (".vtu" in file):
                vtufiles.append(file)
                if "_t_" in file:
                    substring = file.split("_t_")[-1].split(".vtu")[0]
                    if "." in substring:
                        timesteps.append(float(substring.split("_")[0]))
                    else:
                        try:
                            timesteps.append(float(substring.split("_")[0] +
                                               "." +
                                               substring.split("_")[1]))
                        except IndexError:
                            print("""Time step information does not contain
                                          an appropriate separator.""")
                else:
                    print(f"No time information found in vtu file {file}")
        datadict = {"timesteps": timesteps, "vtufiles": vtufiles}
        df = pd.DataFrame(data=datadict)
        df_sorted = df.sort_values(by=["timesteps"], ignore_index=True)
        root = ET.Element("VTKFile")
        root.attrib["type"] = "Collection"
        root.attrib["version"] = "0.1"
        root.attrib["byte_order"] = "LittleEndian"
        root.attrib["compressor"] = "vtkZLibDataCompressor"
        collection = ET.SubElement(root,"Collection")
        timestepselements = []
        #pvdwriter
        for i in df_sorted.index:
            timestepselements.append(ET.SubElement(collection, "DataSet"))
            timestepselements[-1].attrib["timestep"] = str(df_sorted.iloc[i]["timesteps"])
            timestepselements[-1].attrib["group"] = ""
            timestepselements[-1].attrib["part"] = "0"
            timestepselements[-1].attrib["file"] = df_sorted.iloc[i]["vtufiles"]
        tree = ET.ElementTree(root)
        tree.write(filename, encoding="ISO-8859-1",
                   xml_declaration=True, pretty_print=True)

    def write_prj(self, filename):
        """
        exports input data (if available)
        as OGS project files

        Parameters
        ----------
        filename : `str`
        """
        comments = self.tree.xpath("./comment()")
        prjstring = "<OpenGeoSysProject>"
        prjstring += str(comments[0]).split("<OpenGeoSysProject>")[1].split("</OpenGeoSysProject>")[0]
        prjstring += "</OpenGeoSysProject>"
        prjstring = prjstring.replace("\\n","")
        parser = ET.XMLParser(remove_blank_text=True, remove_comments=True)
        root = ET.fromstring(prjstring, parser)
        prjtree = ET.ElementTree(root)
        ET.indent(prjtree, space="    ")
        prjtree.write(filename, encoding="ISO-8859-1", xml_declaration=True, pretty_print=True)

    def write_xdmf(self, filename):
        """
        exports data as XDMF/HDF

        Parameters
        ----------
        filename : `str`
        """
        import meshio
        print("Danger: This function only writes point and cell data. Information could go lost!.")
        mesh = meshio.read(self.vtufilenames[0])
        with meshio.xdmf.TimeSeriesWriter(filename) as writer:
            for i, t in enumerate(self.timesteps):
                mesh = meshio.read(self.vtufilenames[i])
                if i == 0:
                    writer.write_points_cells(mesh.points, mesh.cells)
                writer.write_data(t, point_data=mesh.point_data, cell_data=mesh.cell_data)
