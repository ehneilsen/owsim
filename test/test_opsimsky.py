"""tests for apsupp"""
import unittest

import opsimsky

class TestLSSTSkyModel(unittest.TestCase):
    def setUp(self):
        self.sky_model = opsimsky.SkyModelPre.SkyModelPre(opsimFields=True)

        # Use the first exposures in the baseline2018a sim as a reference
        self.test_mjd = 59853.0167939815
        self.test_mag = 19.3527693482849
        self.test_airmass = 1.00478859066525
        self.test_band = 'z'
        # The fieldID is 1 higher than the field index. Ugh.
        self.test_indx = 1545-1
        
    def test_returnMags(self):
        mag = self.sky_model.returnMags(self.test_mjd, [self.test_indx])
        self.assertAlmostEqual(mag[self.test_band][0], self.test_mag)

    def test_airmass(self):
        airmass = self.sky_model.returnAirmass(self.test_mjd, indx=[self.test_indx])
        self.assertAlmostEqual(airmass[0], self.test_airmass)
        
if __name__ == '__main__':
    unittest.main()
