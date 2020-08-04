import numpy as np
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
        self.points = vtk_to_numpy(output.GetPoints().GetData())
        field = vtk_to_numpy(output.GetPointData().GetArray(fieldname))
        return field
    
    def getFieldnames(self):
        pass

    def readTimeSeries(self,pts = {'pt0': (0.0,0.0,0.0)}):



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
