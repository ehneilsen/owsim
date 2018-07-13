import os

import numpy as np
import astropy.time
import healpy

import lsst.sims.ocs.kernel.time_handler
import lsst.sims.speedObservatory
import lsst.sims.seeingModel
import lsst.sims.cloudModel
import lsst.sims.skybrightness_pre
import lsst.sims.utils 

# constants

SKY_DATA_PATH = '/home/opsim/repos/sims_skybrightness_pre/data'
DEFAULT_START_TIME = astropy.time.Time('2022-01-01 00:00:00Z', format='iso')

# exception classes

# interface functions

calc_five_sigma_depth = lsst.sims.utils.m5_flat_sed
                                
# classes

class SeeingSource(object):

    def __init__(self, survey_start_time=DEFAULT_START_TIME):
        # Separate creation of SeeingSim instance so we do not
        # have to reinstantiate it every iteration
        time_handler = lsst.sims.ocs.kernel.time_handler.TimeHandler(
            survey_start_time.iso.split()[0])
        self.seeing_sim = lsst.sims.seeingModel.SeeingSim(time_handler)

    def __call__(self, seconds_into_survey, band, airmass):
        fwhm500, fwhmGeom, fwhmEff = self.seeing_sim.get_seeing_singlefilter(
            seconds_into_survey, band, airmass)
        return fwhm500, fwhmGeom, fwhmEff
        

class SlewTimeSource(object):
    def __init__(self, lsst_lax_dome=False):
        # Separate creation of Telescope instance so we do not
        # have to reinstantiate it every iteration

        self.telescope = lsst.sims.speedObservatory.Telescope()
        self.lsst_lax_dome = lsst_lax_dome

    def __call__(self, pre_alt, pre_az, pre_band, post_alt, post_az, post_band):
        t = self.telescope.calcSlewTime(
            np.radians(pre_alt), np.radians(pre_az), pre_band,
            np.radians(post_alt), np.radians(post_az), post_band,
            laxDome=self.lsst_lax_dome)
        return t

    
class CloudSource(object):
    def __init__(self, survey_start_time=DEFAULT_START_TIME):
        # Separate creation of CloudModel instance so we do not
        # have to reinstantiate it every iteration
        time_handler = lsst.sims.ocs.kernel.time_handler.TimeHandler(
            survey_start_time.iso.split()[0])
        self.cloud_sim = lsst.sims.cloudModel.CloudModel(time_handler)
        self.cloud_sim.read_data()
        
    def __call__(self, seconds_into_survey):
        cloud_level = self.cloud_sim.get_cloud(seconds_into_survey)
        return cloud_level


class SkyBrightnessSource(object):
    def __init__(self, data_path=SKY_DATA_PATH):
        self.data_path = data_path
        self._sky_sim_hpx = None
        self._sky_sim_field = None
        self._nside = None

    @property
    def sky_sim_hpx(self):
        os.environ['SIMS_SKYBRIGHTNESS_DATA'] = self.data_path
        if self._sky_sim_hpx is None:
            self._sky_sim_hpx = lsst.sims.skybrightness_pre.SkyModelPre(
                opsimFields=False)
        return self._sky_sim_hpx
        
    @property
    def sky_sim_field(self):
        os.environ['SIMS_SKYBRIGHTNESS_DATA'] = self.data_path
        if self._sky_sim_field is None:
            self._sky_sim_field = lsst.sims.skybrightness_pre.SkyModelPre(
                opsimFields=True)
        return self._sky_sim_field
        
    @property
    def nside(self):
        if self._nside is None:
            full = self.sky_sim_hpx.returnMags(
                59853.02, airmass_mask=False, planet_mask=False,
                moon_mask=False, zenith_mask=False)
            self._nside = healpy.npix2nside(full["r"].size)
        return self._nside
        
    def field_brightness(self, mjd, idx, band):
        # sky indexes are 0 indexed, field indexes are 1 indexed
        sky_idx = idx-1
        sky_mag = self.sky_sim_field.returnMags(mjd, [sky_idx])[band]
        return sky_mag

    def hpx_brightness(self, mjd, idx, band):
        sky_mag = self.sky_sim_hpx.returnMags(mjd, [idx])[band]
        return sky_mag

    def coord_brightness(self, mjd, ra, decl, band):
        theta = np.radians(90-decl)
        phi = np.radians(ra)
        hpix = healpy.ang2pix(self.nside, theta, phi)
        sky_mag = self.hpx_brightness(mjd, hpix, band)
        return sky_mag

# internal functions & classes
