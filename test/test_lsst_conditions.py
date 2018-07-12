"""tests for lsst_conditions"""

import numpy as np
import unittest
from lsst_conditions import SeeingSource
from lsst_conditions import SlewTimeSource
from lsst_conditions import CloudSource
from lsst_conditions import SkyBrightnessSource

class TestSeeingSource(unittest.TestCase):
    def test_single(self):
        seeing_source = SeeingSource()
        fwhm_500, fwhm_geom, fwhm_eff = seeing_source(42, 'g', 1.2)
        self.assertTrue(0.2 < fwhm_500 < 2.5)
        self.assertTrue(0.2 < fwhm_geom < 2.5)
        self.assertTrue(0.2 < fwhm_eff < 2.5)
        
    def test_vectorized(self):
        seeing_source = np.vectorize(SeeingSource())
        idx = np.array([2, 4303, 35603, 44221])
        band = np.array(['g', 'r', 'g', 'z'])
        airmass = np.array([1.02, 2.2, 1.3, 1.2])
                       
        fwhm_500, fwhm_geom, fwhm_eff = seeing_source(idx, band, airmass)
        for i in range(len(idx)):
            self.assertTrue(0.2 < fwhm_500[i] < 2.5)
            self.assertTrue(0.2 < fwhm_geom[i] < 2.5)
            self.assertTrue(0.2 < fwhm_eff[i] < 2.5)

class TestSlewTimeSource(unittest.TestCase):
    def test_single(self):
        slew_time_source = SlewTimeSource()
        t = slew_time_source(80.1, 44.2, 'g',
                             75.2, 50.1, 'z')
        self.assertTrue(0 < t < 300)

    def test_vectorized(self):
        npts = 5
        bands = ['u', 'g', 'r', 'i', 'z']
        slew_time_source = np.vectorize(SlewTimeSource())
        pre_alt = 30 + np.random.rand(npts) * 59
        pre_az = np.random.rand(npts) * 360
        pre_filter = np.random.choice(bands, npts)
        post_alt = 30 + np.random.rand(npts) * 59
        post_az = np.random.rand(npts) * 360
        post_filter = np.random.choice(bands, npts)
        t = slew_time_source(pre_alt, pre_az, pre_filter,
                             post_alt, post_az, post_filter)
        self.assertEqual(len(t), npts)
        self.assertTrue(np.all(0 < t))
        self.assertTrue(np.all(t < 300))


class TestCloudSource(unittest.TestCase):
    def test_single(self):
        cloud_source = CloudSource()
        clouds = cloud_source(42.234 * 24*3600)
        self.assertTrue(0 <= clouds < 9)

    def test_vectorized(self):
        npts = 5
        cloud_source = np.vectorize(CloudSource())
        seconds_into_survey = np.random.rand(npts)*(24*60*60*265*10)
        clouds = cloud_source(seconds_into_survey)
        self.assertEqual(len(clouds), npts)
        self.assertTrue(np.all(0 <= clouds))
        self.assertTrue(np.all(clouds < 9))


class TestSkyBrightnessSource(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestSkyBrightnessSource, self).__init__(*args, **kwargs)
        self.sky_source = SkyBrightnessSource()
    
    def test_single_field_brightness(self):
        sky = self.sky_source.field_brightness(
            59853.0167939815, 1545, 'z')
        self.assertAlmostEqual(sky, 19.3528, delta=0.05)
                                          
    def test_single_hpx_brightness(self):
        sky = self.sky_source.hpx_brightness(
            59853.0167939815, 8748, 'z')
        self.assertAlmostEqual(sky, 19.3528, delta=0.05)

    def test_single_coord_brightness(self):
        sky = self.sky_source.coord_brightness(
            59853.0167939815, 305.088793, -24.889283, 'z')
        self.assertAlmostEqual(sky, 19.3528, delta=0.05)

        

    
if __name__ == '__main__':
    unittest.main()
