import unittest
import numpy as np
from context import vtuIO

class TestiOGS(unittest.TestCase):

    def test_pvd_read(self):
        pvdfile=vtuIO.PVDIO("../examples", "square_1e2_pcs_0.pvd", dim=2)
        time = pvdfile.timesteps
        selected_points = {'pt0': (0.25, 0.5, 0.0), 'pt1': (0.75, 0.5, 0.0)}
        pressure_interpolation = pvdfile.read_time_series('pressure', selected_points)
        self.assertAlmostEqual(pressure_interpolation['pt0'][-1],0.5)
        self.assertAlmostEqual(pressure_interpolation['pt1'][-1],-0.5)

    def test_point_set_read(self):
        t = 0.5
        pvdfile=vtuIO.PVDIO("../examples", "square_1e2_pcs_0.pvd", dim=2)
        xaxis =  [(i,0,0) for i in np.linspace(start=0.0, stop=1.0, num=100)]
        y_pred = np.linspace(start=0.5, stop=-0.5, num=100)
        pressure_xaxis_t1 = pvdfile.read_point_set_data(t, 'pressure', pointsetarray=xaxis)
        for i, p in enumerate(pressure_xaxis_t1):
            self.assertAlmostEqual(y_pred[i],p)

    def test_vtu_write_read(self):
        vtufile = vtuIO.VTUIO("../examples/square_1e2_pcs_0_ts_1_t_1.000000.vtu", dim=2)
        def fct(x,y,z):
            return x*10
        def fct2(x,y,z):
            return -y*10
        vtufile.func_to_field(fct, "field1","fields.vtu")
        vtufile = vtuIO.VTUIO("fields.vtu", dim=2)
        vtufile.func_to_m_dim_field([fct,fct2], "field2","fields.vtu")
        vtufile = vtuIO.VTUIO("fields.vtu", dim=2)
        f1 = vtufile.get_point_data("field1", pts={'pt0':(0.75,0.0,0.0)})
        self.assertAlmostEqual(f1['pt0'],7.5)
        f2 = vtufile.get_point_data("field2", pts={'pt0':(0.25,0.25,0.0)})
        self.assertAlmostEqual(f2['pt0'][1], -2.5)
        self.assertAlmostEqual(f2['pt0'][0], 2.5)


if __name__ == '__main__':
    unittest.main()
