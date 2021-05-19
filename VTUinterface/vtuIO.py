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
import numpy as np
import pandas as pd
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
    """
    def __init__(self, filename, nneighbors=20, dim=3):
        self.filename = filename
        self.reader = vtk.vtkXMLUnstructuredGridReader()
        self.reader.SetFileName(self.filename)
        self.reader.Update()
        self.output = self.reader.GetOutput()
        self.pdata = self.output.GetPointData()
        self.points = vtk_to_numpy(self.output.GetPoints().GetData())
        self.dim = dim
        self.nneighbors = nneighbors
        if self.dim == 2:
            self.points = np.delete(self.points, 2, 1)

    def getNeighbors(self, points_interpol):
        """
        Method for obtaining neighbor points for interpolation.
        """
        df = pd.DataFrame(self.points)
        neighbors = {}
        for i, (_, val) in enumerate(points_interpol.items()):
            if self.dim == 2:
                df["r_"+str(i)] = (df[0]-val[0]) * (df[0]-val[0]) + (df[1]-val[1]) * (df[1]-val[1])
            else:
                df["r_"+str(i)] = ((df[0]-val[0]) * (df[0]-val[0]) + (df[1]-val[1]) * (df[1]-val[1])
                        + (df[2]-val[2]) * (df[2]-val[2]))
            neighbors[i] = df.sort_values(by=["r_" + str(i)]).head(self.nneighbors).index
        return neighbors

    def getData(self, neighbors, points_interpol, fieldname, interpolation_method="linear"):
        """
        Get interpolated data for points_interpol using neighbor points.
        """
        field = self.getField(fieldname)
        resp = {}
        for i, (key, val) in enumerate(points_interpol.items()):
            if self.dim == 1:
                data = pd.DataFrame(self.points[:,0], columns = ['x'])
                data["y"] = field
                data.sort_values(by = ['x'], inplace=True)
                data.drop_duplicates(subset=['x'], inplace=True)
                f = interp1d(data['x'], data['y'])
                resp[key] = f(val[0])
            elif self.dim == 2:
                grid_x, grid_y = np.array([[[val[0]]],[[val[1]]]])
                resp[key] = griddata(self.points[neighbors[i]], field[neighbors[i]],
                        (grid_x, grid_y), method=interpolation_method)[0][0]
            else:
                grid_x, grid_y, grid_z = np.array([[[[val[0]]]], [[[val[1]]]], [[[val[2]]]]])
                resp[key] = griddata(self.points[neighbors[i]], field[neighbors[i]],
                        (grid_x, grid_y, grid_z), method=interpolation_method)[0][0][0]
        return resp

    def getField(self, fieldname):
        """
        Return vtu point field as numpy array.
        fieldname : `str`
        """
        field = vtk_to_numpy(self.pdata.GetArray(fieldname))
        return field

    def getFieldNames(self):
        """
        Get names of all point fields in the vtu file.
        """
        fieldnames = []
        for i in range(self.pdata.GetNumberOfArrays()):
            fieldnames.append(self.pdata.GetArrayName(i))
        return fieldnames

    def getPointData(self, fieldname, pts = None, interpolation_method="linear"):
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
        nb = self.getNeighbors(pts)
        if isinstance(fieldname, str):
            data = self.getData(nb, pts, fieldname, interpolation_method=interpolation_method)
            for pt in pts:
                resp[pt]=data[pt]
        elif isinstance(fieldname, list):
            data = {}
            for field in fieldname:
                data[field] = self.getData(nb, pts, field, interpolation_method=interpolation_method)
            for pt in pts:
                for field in fieldname:
                    resp[pt][field] = data[field][pt]
        return resp

    def getPointSetData(self, fieldname, pointsetarray=None, interpolation_method="linear"):
        """
        Get data specified in fieldname at all points specified in "pointsetarray".

        Parameters
        ----------
        fieldname : `str`
        pointsetarray : `list`, optional
                        default: [(0,0,0)]
        interpolation_method : `str`, optional
                               default: 'linear'
        """
        if pointsetarray is None:
            pointsetarray = [(0,0,0)]
        pts = {}
        # convert into point dictionary
        for i, entry in enumerate(pointsetarray):
            pts['pt'+str(i)] = entry
        resp = self.getPointData(fieldname, pts=pts, interpolation_method=interpolation_method)
        resp_list = []
        # convert point dictionary into list
        for i, entry in enumerate(pointsetarray):
            resp_list.append(resp['pt' + str(i)])
        resp_array = np.array(resp_list)
        return resp_array

    def func2Field(self, function, fieldname, ofilename, cell=False):
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
        fieldarray = np.zeros(len(self.points))
        for i,_ in enumerate(fieldarray):
            if self.dim == 2:
                fieldarray[i] = function(self.points[i,0], self.points[i,1], 0.0)
            else:
                fieldarray[i] = function(self.points[i,0], self.points[i,1], self.points[i,2])
        field_vtk = numpy_to_vtk(fieldarray)
        r = self.pdata.AddArray(field_vtk)
        self.pdata.GetArray(r).SetName(fieldname)
        if cell is True:
            p2c = vtk.vtkPointDataToCellData()
            p2c.SetInputData(self.output)
            p2c.Update()
            outcells = p2c.GetOutput()
            cells = outcells.GetCellData()
            array =  cells.GetArray(fieldname)
            cells_orig = self.output.GetCellData()
            cells_orig.AddArray(array)
            self.pdata.RemoveArray(fieldname)
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(ofilename)
        writer.SetInputData(self.output)
        writer.Write()

    def func2MdimField(self, functionarray, fieldname, ofilename, cell=False):
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
        fieldarray = np.zeros((len(self.points), mdim))
        for i,_ in enumerate(fieldarray):
            for j, func in enumerate(functionarray):
                if self.dim == 2:
                    fieldarray[i,j] = func(
                        self.points[i,0],
                        self.points[i,1],
                        0.0)
                else:
                    fieldarray[i,j] = func(
                        self.points[i,0],
                        self.points[i,1],
                        self.points[i,2])
        field_vtk = numpy_to_vtk(fieldarray)
        r = self.pdata.AddArray(field_vtk)
        self.pdata.GetArray(r).SetName(fieldname)
        if cell is True:
            p2c = vtk.vtkPointDataToCellData()
            p2c.SetInputData(self.output)
            p2c.Update()
            outcells = p2c.GetOutput()
            cells = outcells.GetCellData()
            array =  cells.GetArray(fieldname)
            cells_orig = self.output.GetCellData()
            cells_orig.AddArray(array)
            self.pdata.RemoveArray(fieldname)
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(ofilename)
        writer.SetInputData(self.output)
        writer.Write()

    def pointData2CellData(self, fieldname, ofilename):
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
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(ofilename)
        writer.SetInputData(self.output)
        writer.Write()

    def writeField(self, field, fieldname, ofilename):
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
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(ofilename)
        writer.SetInputData(self.output)
        writer.Write()



class PVDIO:
    """
    Class for handling I/O of PVD files

    Parameters
    ----------
    folder : `str`
    filename : `str`
    nneighbors : `int`, optional
    dim : `int`
    """
    def __init__(self, folder, filename, nneighbors=20, dim=3):
        self.folder = folder
        self.filename = ""
        self.nneighbors = nneighbors
        self.timesteps = np.array([])
        self.vtufilenames = []
        self.readPVD(os.path.join(folder, filename))
        self.dim = dim

    def readPVD(self, filename):
        """
        Read in PVD file

        Parameters
        ----------
        filename : `str`
        """
        print(filename)
        self.filename = filename
        tree = ET.parse(self.filename)
        root = tree.getroot()
        for collection in root.getchildren():
            for dataset in collection.getchildren():
                self.timesteps = np.append(self.timesteps, [float(dataset.attrib['timestep'])])
                self.vtufilenames.append(dataset.attrib['file'])

    def readTimeSeries(self, fieldname, pts=None, interpolation_method="linear"):
        """
        Return time series data of field "fieldname" at points pts.
        Also a list of fieldnames can be provided as "fieldname"

        Parameters
        ----------
        fieldname : `str`
        pts : `dict`, optional
              default: {'pt0': (0.0,0.0,0.0)}
        interpolation_method : `str`, optional
                               default: 'linear
        """
        if pts is None:
            pts = {'pt0': (0.0,0.0,0.0)}
        resp_t = {}
        for pt in pts:
            if isinstance(fieldname, str):
                resp_t[pt] = []
            elif isinstance(fieldname, list):
                resp_t[pt] = {}
                for field in fieldname:
                    resp_t[pt][field] = []
        for i, filename in enumerate(self.vtufilenames):
            # quick and dirty trick for serial pvtu files:
            # TODO: real handling of parallel files
            fn_new = filename.replace(".pvtu", "_0.vtu")
            vtu = VTUIO(os.path.join(self.folder,fn_new),
                    nneighbors=self.nneighbors, dim=self.dim)
            if i == 0:
                nb = vtu.getNeighbors(pts)
            if isinstance(fieldname, str):
                data = vtu.getData(nb, pts, fieldname, interpolation_method=interpolation_method)
                for pt in pts:
                    resp_t[pt].append(data[pt])
            elif isinstance(fieldname, list):
                data = {}
                for field in fieldname:
                    data[field] = vtu.getData(nb, pts, field, interpolation_method=interpolation_method)
                for pt in pts:
                    for field in fieldname:
                        resp_t[pt][field].append(data[field][pt])
        resp_t_array = {}
        for pt, field in resp_t.items():
            if isinstance(fieldname, str):
                resp_t_array[pt] = np.array(field)
            elif isinstance(fieldname, list):
                resp_t_array[pt] = {}
                for field_, fieldarray in resp_t[pt].items():
                    resp_t_array[pt][field_] = np.array(fieldarray)
        return resp_t_array

    def readTimeStep(self, timestep, fieldname):
        """
        Print field "fieldname" at time "timestep".

        Parameters
        ----------
        timestep : `int`
        fieldname : `str`
        """
        filename = None
        for i, ts in enumerate(self.timesteps):
            if timestep == ts:
                filename = self.vtufilenames[i]
        if not filename is None:
            vtu = VTUIO(os.path.join(self.folder,filename),
                    nneighbors=self.nneighbors, dim=self.dim)
            field = vtu.getField(fieldname)
        else:
            filename1 = None
            filename2 = None
            timestep1 = 0.0
            timestep2 = 0.0
            for i, ts in enumerate(self.timesteps):
                try:
                    if (timestep > ts) and (timestep < self.timesteps[i+1]):
                        timestep1 = ts
                        timestep2 = self.timesteps[i+1]
                        filename1 = self.vtufilenames[i]
                        filename2 = self.vtufilenames[i+1]
                except IndexError:
                    print("time step is out of range")
            if (filename1 is None) or (filename2 is None):
                print("time step is out of range")
            else:
                vtu1 = VTUIO(os.path.join(self.folder,filename1),
                        nneighbors=self.nneighbors, dim=self.dim)
                vtu2 = VTUIO(os.path.join(self.folder,filename2),
                        nneighbors=self.nneighbors, dim=self.dim)
                field1 = vtu1.getField(fieldname)
                field2 = vtu2.getField(fieldname)
                fieldslope = (field2-field1)/(timestep2-timestep1)
                field = field1 + fieldslope * (timestep-timestep1)
        return field

    def readPointSetData(self, timestep, fieldname, pointsetarray = None, interpolation_method="linear"):
        """
        Get data of field "fieldname" at time "timestep" alon a given "pointsetarray".

        Parameters
        ----------
        timestep : `int`
        fieldname : `str`
        pointsetarray : `array`, optional
                        default: [(0,0,0)]
        interpolation_method : `str`
                               default: 'linear'
        """
        if pointsetarray is None:
            pointsetarray = [(0,0,0)]
        filename = None
        for i, ts in enumerate(self.timesteps):
            if timestep == ts:
                filename = self.vtufilenames[i]
        if not filename is None:
            vtu = VTUIO(os.path.join(self.folder,filename),
                    nneighbors=self.nneighbors, dim=self.dim)
            field = vtu.getPointSetData(fieldname, pointsetarray, interpolation_method=interpolation_method)
        else:
            filename1 = None
            filename2 = None
            timestep1 = 0.0
            timestep2 = 0.0
            for i, ts in enumerate(self.timesteps):
                try:
                    if (timestep > ts) and (timestep < self.timesteps[i+1]):
                        timestep1 = ts
                        timestep2 = self.timesteps[i+1]
                        filename1 = self.vtufilenames[i]
                        filename2 = self.vtufilenames[i+1]
                except IndexError:
                    print("time step is out of range")
            if (filename1 is None) or (filename2 is None):
                print("time step is out of range")
            else:
                vtu1 = VTUIO(os.path.join(self.folder,filename1),
                    nneighbors=self.nneighbors, dim=self.dim)
                vtu2 = VTUIO(os.path.join(self.folder,filename2),
                    nneighbors=self.nneighbors, dim=self.dim)
                field1 = vtu1.getPointSetData(fieldname, pointsetarray, interpolation_method=interpolation_method)
                field2 = vtu2.getPointSetData(fieldname, pointsetarray, interpolation_method=interpolation_method)
                fieldslope = (field2-field1)/(timestep2-timestep1)
                field = field1 + fieldslope * (timestep-timestep1)
        return field

    def clearPVDRelPath(self):
        """
        Delete relative directory paths in the vtu filenames of the PVD file.
        """
        xpath="./Collection/DataSet"
        tree = ET.parse(self.filename)
        root = tree.getroot()
        find_xpath = root.findall(xpath)
        for tag in find_xpath:
            filename = tag.get("file")
            filename_new = filename.split("/")[-1]
            tag.set("file", filename_new)
        tree.write(self.filename,
                            encoding="ISO-8859-1",
                            xml_declaration=True,
                            pretty_print=True)
        #update file list:
        newlist = []
        for entry in  self.vtufilenames:
            newlist.append(entry.split("/")[-1])
        self.vtufilenames = newlist
