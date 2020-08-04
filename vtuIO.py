import numpy as np
import pandas as pd
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
from lxml import etree as ET
from scipy.interpolate import griddata


class VTUIO(object):
    def __init__(self):
        self.filename = ""
        self.ts_files = {}
        self.ts_files['ts'] = []
        self.ts_files['filename'] = []
        self.points = np.array([])

    def _getNeighbors(self, points_interpol, df):
        neighbors = {}
        for i, key in enumerate(points_interpol):
            df["r_"+str(i)]=(df[0]-points_interpol[key][0])*(df[0]-points_interpol[key][0])+(df[1]-points_interpol[key][1])*(df[1]-points_interpol[key][1])+(df[2]-points_interpol[key][2])*(df[key]-points_interpol[key][2])
            neighbors[i] = df.sort_values(by=["r_"+str(i)]).head(10).index
        return neighbors
    def _getData(self, neighbors, pts, filename, fieldname):
        field = self.readVTUfile(filename,fieldname)
        resp = {}
        for i, key in enumerate(points_interpol):
            grid_x, grid_y, grid_z = np.mgrid[points_interpol[key][0]:(points_interpol[key][0]+0.1):1, points_interpol[key][1]:(points_interpol[key][1]+0.1):1, points_interpol[key][2]:(points_interpol[key][2]+0.1):]
            resp[key] = griddata(self.points[neighbors[i]], field[neighbors[i]], (grid_x, grid_y, grid_z), method='linear')[0][0][0]
        return resp

    def readPVD(self,filename, prefix='./'):
        self.filename = prefix + filename
        tree = ET.parse(self.filename)
        root = tree.getroot()
        for collection in root.getchildren():
            for dataset in collection.getchildren():
                self.ts_files['ts'].append(dataset.attrib['timestep'])
                self.ts_files['filename'].append(dataset.attrib['file'])

    def readTimeStep(self, timestep, fieldname):
        if not len(self.ts_files['ts']) > 0:
            raise LookupError
        else:
            for i in self.ts_files['ts']:
                if timestep == i:
                    filename = self.ts_files['filename'][i]
        field = self.readVTUfile(filename, fieldname)
        return field

    def readVTUfile(self, filename, fieldname):
        reader = vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()
        output = reader.GetOutput()
        if not len(self.points) > 0
            self.points = vtk_to_numpy(output.GetPoints().GetData())
        field = vtk_to_numpy(output.GetPointData().GetArray(fieldname))
        return field
    
    def getFieldnames(self):
        pass

    def readTimeSeries(self,pts = {'pt0': (0.0,0.0,0.0)},fieldname):
        df = pd.DataFrame(self.points)
        nb = self._getNeighbors(pts, df)
        resp_t = {}
        for pt in pts:
            resp_t[pt] = []
        for filename in self.ts_files['filename']
            data = _getData(nb, pts, filename, fieldname)
        for pt in pts:
            resp_t[pt].append(data[pt])
        return resp_t

    def readVTU(self,filename,observation_points):
        reader = vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()
        output = reader.GetOutput()
        points = vtk_to_numpy(output.GetPoints().GetData())
        delta=np.zeros(len(points))
        indices = {}
        response = {}
        temp = vtk_to_numpy(output.GetPointData().GetArray('temperature_interpolated'))
        press = vtk_to_numpy(output.GetPointData().GetArray('pressure_interpolated'))
        displ = vtk_to_numpy(output.GetPointData().GetArray('displacement'))
        stress = vtk_to_numpy(output.GetPointData().GetArray('sigma'))
        for point in observation_points:
            delta[:] = ((observation_points[point][0]-points[:,0])**2
                    +(observation_points[point][1]-points[:,1])**2
                    +(observation_points[point][2]-points[:,2])**2)
            indices[point] = delta.argmin()
            response[point] = {
                    'x_grid': points[indices[point],0],
                    'y_grid': points[indices[point],1],
                    'z_grid': points[indices[point],2],
                    'temp': temp[indices[point]],
                    'press': press[indices[point]],
                    'ux': displ[indices[point],0],
                    'uy': displ[indices[point],1],
                    'sigmaxx': stress[indices[point],0],
                    'sigmayy': stress[indices[point],1]}
~~~~~~~~~~~~
        return response
