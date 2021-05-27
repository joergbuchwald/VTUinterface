import unittest
import numpy as np
from context import VTUinterface

class TestiOGS(unittest.TestCase):

    def test_pvd_read(self):
        pvdfile = VTUinterface.PVDIO("examples/square_1e2_pcs_0.pvd", dim=2)
        time = pvdfile.timesteps
        selected_points = {'pt0': (0.25, 0.5, 0.0), 'pt1': (0.75, 0.5, 0.0)}
        pressure_interpolation = pvdfile.read_time_series('pressure', selected_points)
        self.assertAlmostEqual(pressure_interpolation['pt0'][-1],0.5)
        self.assertAlmostEqual(pressure_interpolation['pt1'][-1],-0.5)

    def test_point_set_read(self):
        t = 0.5
        pvdfile = VTUinterface.PVDIO("examples/square_1e2_pcs_0.pvd", dim=2)
        xaxis =  [(i,0,0) for i in np.linspace(start=0.0, stop=1.0, num=100)]
        y_pred = np.linspace(start=0.5, stop=-0.5, num=100)
        pressure_xaxis_t1 = pvdfile.read_point_set_data(t, 'pressure', pointsetarray=xaxis)
        for i, p in enumerate(pressure_xaxis_t1):
            self.assertAlmostEqual(y_pred[i],p)

    def test_vtu_write_read(self):
        vtufile = VTUinterface.VTUIO("examples/square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)
        def fct(x,y,z):
            return x*10
        def fct2(x,y,z):
            return -y*10
        vtufile.func_to_field(fct, "field1","fields.vtu")
        vtufile = VTUinterface.VTUIO("fields.vtu", dim=2)
        vtufile.func_to_m_dim_field([fct,fct2], "field2","fields.vtu")
        vtufile = VTUinterface.VTUIO("fields.vtu", dim=2)
        f1 = vtufile.get_point_data("field1", pts={'pt0':(0.75,0.0,0.0)})
        self.assertAlmostEqual(f1['pt0'],7.5)
        f2 = vtufile.get_point_data("field2", pts={'pt0':(0.25,0.25,0.0)})
        self.assertAlmostEqual(f2['pt0'][1], -2.5)
        self.assertAlmostEqual(f2['pt0'][0], 2.5)

    def test_read_1d_field(self):
        vtufile = VTUinterface.VTUIO("examples/line_1_time_dep_dirichlet.vtu", dim=1)
        vtupflist = vtufile.get_point_field_names()
        pflist = [f"t_{i+1}s" for i in range(10)]
        for i, entry in enumerate(vtupflist):
            self.assertEqual(entry, pflist[i])
        pts = {'pt0': (0.33,0,0), 'pt1': (0.97,0,0)}
        data = vtufile.get_point_data(pflist, pts=pts)
        for pt in pts:
            for i, field in enumerate(pflist):
                self.assertAlmostEqual(float(data[pt][field]),(i+1)*pts[pt][0])
    def test_point_to_celldata(self):
        vtufile = VTUinterface.VTUIO("examples/line_1_time_dep_dirichlet.vtu", dim=1)
        vtufile.point_data_to_cell_data("t_10s", "line_1_time_dep_dirichlet_cdata.vtu")
        vtufile = VTUinterface.VTUIO("line_1_time_dep_dirichlet_cdata.vtu", dim=1)
        cflist = ["t_10s"]
        self.assertEqual(len(cflist),len(vtufile.get_cell_field_names()))
        self.assertEqual(vtufile.get_cell_field_names()[0], cflist[0])
        field = vtufile.get_cell_field_as_point_data("t_10s")
        for i, entry in enumerate(field):
            if i == 0:
                self.assertEqual(vtufile.points[i], 0)
                self.assertEqual(entry, 0.5)
            elif i == 10:
                self.assertEqual(vtufile.points[i], 1)
                self.assertEqual(entry, 9.5)
            else:
                self.assertAlmostEqual(entry, vtufile.points[i]*10)
    def test_read_time_step(self):
        t1 = 0.5
        t2 = 1
        pvdfile = VTUinterface.PVDIO("examples/square_1e2_pcs_0.pvd", dim=2)
        field_last_step = pvdfile.read_time_step(t2, "pressure")
        field_t1 = pvdfile.read_time_step(t1, "pressure")
        for i, entry in enumerate(field_t1):
            self.assertAlmostEqual(0.5*field_last_step[i], entry)

if __name__ == '__main__':
    unittest.main()
