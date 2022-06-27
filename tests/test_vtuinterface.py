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

    def test_pvd_read_vtk(self):
        pvdfile = VTUinterface.PVDIO("examples/square_1e2_pcs_0.pvd", dim=2, interpolation_backend="vtk")
        time = pvdfile.timesteps
        selected_points = {'pt0': (0.25, 0.5, 0.0), 'pt1': (0.75, 0.5, 0.0)}
        pressure_interpolation_v = pvdfile.read_time_series('pressure', selected_points, interpolation_method="voronoi")
        # not equal to test test_pvd_read():
        self.assertAlmostEqual(pressure_interpolation_v['pt0'][-1],0.6)
        self.assertAlmostEqual(pressure_interpolation_v['pt1'][-1],-0.6)
        pressure_interpolation_s = pvdfile.read_time_series('pressure', selected_points, interpolation_method="shepard")
        self.assertAlmostEqual(pressure_interpolation_s['pt0'][-1], 0.41946054)
        self.assertAlmostEqual(pressure_interpolation_s['pt1'][-1],-0.41946054)

    def test_point_set_read(self):
        t = 0.5
        pvdfile = VTUinterface.PVDIO("examples/square_1e2_pcs_0.pvd", dim=2)
        xaxis =  [(i,0,0) for i in np.linspace(start=0.0, stop=1.0, num=100)]
        y_pred = np.linspace(start=0.5, stop=-0.5, num=100)
        pressure_xaxis_t1 = pvdfile.read_set_data(t, 'pressure', data_type="point", pointsetarray=xaxis)
        for i, p in enumerate(pressure_xaxis_t1):
            self.assertAlmostEqual(y_pred[i],p)

    def test_vtu_cell_func_write_read(self):
        vtufile = VTUinterface.VTUIO("examples/square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)
        field = vtufile.get_point_field("pressure")
        fieldnew = 0.5*field
        vtufile.write_point_field(fieldnew, "pressure_new","write_test.vtu")
        vtufile = VTUinterface.VTUIO("write_test.vtu", dim=2)
        def fct(x,y,z):
            return x*10
        def fct2(x,y,z):
            return -y*10
        vtufile.func_to_field(fct, "field1","write_test.vtu", cell=True)
        vtufile = VTUinterface.VTUIO("write_test.vtu", dim=2)
        vtufile.func_to_m_dim_field([fct,fct2], "field2","write_test.vtu", cell=True)
        vtufile = VTUinterface.VTUIO("write_test.vtu", dim=2)
        self.assertTrue("pressure_new" in vtufile.get_point_field_names())
        self.assertTrue("field1" in vtufile.get_cell_field_names())
        self.assertTrue("field2" in vtufile.get_cell_field_names())

    def test_vtu_func_write_read(self):
        vtufile = VTUinterface.VTUIO("examples/square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)
        def fct(x,y,z):
            return x*10
        def fct2(x,y,z):
            return -y*10
        vtufile.func_to_field(fct, "field1","fields.vtu")
        vtufile = VTUinterface.VTUIO("fields.vtu", dim=2)
        vtufile.func_to_m_dim_field([fct,fct2], "field2","fields.vtu")
        vtufile = VTUinterface.VTUIO("fields.vtu", dim=2)
        f1 = vtufile.get_data("field1", pts={'pt0':(0.75,0.0,0.0)})
        self.assertAlmostEqual(f1['pt0'],7.5)
        f2 = vtufile.get_data("field2", pts={'pt0':(0.25,0.25,0.0)})
        self.assertAlmostEqual(f2['pt0'][1], -2.5)
        self.assertAlmostEqual(f2['pt0'][0], 2.5)

    def test_read_1d_field(self):
        vtufile = VTUinterface.VTUIO("examples/line_1_time_dep_dirichlet.vtu", dim=1)
        vtupflist = vtufile.get_point_field_names()
        pflist = [f"t_{i+1}s" for i in range(10)]
        for i, entry in enumerate(vtupflist):
            self.assertEqual(entry, pflist[i])
        pts = {'pt0': (0.33,0,0), 'pt1': (0.97,0,0)}
        data = vtufile.get_data(pflist, pts=pts)
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
    def test_celldata_as_point_data(self):
        vtufile = VTUinterface.VTUIO("examples/square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)
        v = vtufile.get_cell_field_as_point_data("MaterialIDs")
        self.assertEqual(v[0], 0.0)
    def test_center_points(self):
        vtufile = VTUinterface.VTUIO("examples/square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)
        points = vtufile.cell_center_points
        self.assertEqual(points[0].all(), np.array([0.05, 0.05]).all())
    def test_time_series_cell(self):
        pvdfile = VTUinterface.PVDIO("examples/square_1e2_pcs_0.pvd", dim=2)
        ts = pvdfile.read_time_series("MaterialIDs", pts={"pt0":[0.345,0.5231,0]}, data_type="cell")
        self.assertEqual(ts["pt0"].all(), np.array([0.0, 0.0]).all())
    def test_read_time_step(self):
        t1 = 0.5
        t2 = 1
        pvdfile = VTUinterface.PVDIO("examples/square_1e2_pcs_0.pvd", dim=2)
        field_last_step = pvdfile.read_time_slice(t2, "pressure")
        field_t1 = pvdfile.read_time_slice(t1, "pressure")
        for i, entry in enumerate(field_t1):
            self.assertAlmostEqual(0.5*field_last_step[i], entry)
    def test_read_pvtu(self):
        f = VTUinterface.PVDIO("examples/run_0_results.pvd", dim=2)
        f.clear_pvd_rel_path(write=False)
        p = f.read_time_series("pressure", {"pt0":(1.53,1.73,0)})
        self.assertAlmostEqual(p["pt0"].all(), np.array([  100000., 11214944.35401228]).all())

    def test_nearest_points(self):
        vtufile = VTUinterface.VTUIO("examples/square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)
        pts = {"pt0": (0.07,0.07,0.0),"pt1": (0.02,0.02,0.0),"pt2":(0.02,0.07,0.0), "pt3": (0.07,0.02,0.0)}
        points = vtufile.get_nearest_points(pts)
        indices = vtufile.get_nearest_indices(pts)
        self.assertEqual(indices["pt0"], 12)
        self.assertEqual(indices["pt1"], 0)
        self.assertEqual(indices["pt2"], 11)
        self.assertEqual(indices["pt3"], 1)
        self.assertEqual(points["pt0"].all(), np.array([0.1, 0.1]).all())
        self.assertEqual(points["pt1"].all(), np.array([0.0, 0.0]).all())
        self.assertEqual(points["pt2"].all(), np.array([0.0, 0.1]).all())
        self.assertEqual(points["pt3"].all(), np.array([0.1, 0.0]).all())
    def test_aggregate(self):
        pvdfile1 = VTUinterface.PVDIO("examples/tunnel_heat_tunnel_restart.pvd", dim=2)
        pvdfile2 = VTUinterface.PVDIO("examples/tunnel_heat_tunnel_inner.pvd", dim=2)
        T_max_heater1 = pvdfile1.read_aggregate("temperature",agg_fct="max",
                pointsetarray="examples/tunnel_heat_tunnel_inner_ts_160_t_9856003.000000.vtu")
        T_max_heater2 = pvdfile2.read_aggregate("temperature",agg_fct="max")
        self.assertEqual(np.array(T_max_heater1).all(), np.array(T_max_heater2).all())

if __name__ == '__main__':
    unittest.main()
